from __future__ import annotations

import random
import sys
from datetime import datetime
from pathlib import Path

import dask
import numpy as np
import tifffile
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
# PyTorch TensorBoard support
from torch.utils.tensorboard import SummaryWriter
from torchvision.transforms import transforms
from tqdm import tqdm

dask.config.set(scheduler='synchronous')


def partial_cellpose_loss(lbl, y, device, mask_gradient):
    """
    Calculates the loss function between true labels lbl and prediction y.
    Args:
        lbl (numpy.ndarray): True labels (cellprob, flowsY, flowsX).
        y (torch.Tensor): Predicted values (flowsY, flowsX, cellprob).
        device (torch.device): Device on which the tensors are located.
    Returns:
        torch.Tensor: Loss value.
    """

    criterion = nn.MSELoss(reduction="none")  # Change to none to keep all losses
    criterion2 = nn.BCEWithLogitsLoss(reduction="none")  # Change to none to keep all losses

    assert len(lbl[0]) == 3
    assert len(y[0]) == 3
    if lbl.device != device:
        lbl = lbl.to(device)
    veci = 5. * lbl[:, 1:]

    loss = criterion(y[:, :2], veci)
    loss_y = loss[:, 0][mask_gradient > 0]
    loss_x = loss[:, 1][mask_gradient > 0]
    loss = loss_y.mean() + loss_x.mean()
    loss /= 2.

    loss2 = criterion2(y[:, -1], (lbl[:, 0] > 0.5).float())  # not optimal
    loss2 = loss2[mask_gradient > 0]
    loss2 = loss2.mean()
    loss = loss + loss2

    return loss


def train_one_epoch(
        device: torch.device | str,
        epoch_index: int,
        rna2seg: nn.Module,
        training_loader: torch.utils.data.dataloader.DataLoader,
        optimizer: torch.optim.Optimizer,
        print_loss_every=4,
        tb_writer: SummaryWriter = None,
        validation_loader: torch.utils.data.dataloader.DataLoader = None,
        path_save_model: str | Path = None,
        cellbound_prob: float = 0.8,
        no_cellbound_first_iter: int = 0,
        cellbound_index_to_use: list[int] | None = None,
        best_val_loss: float | None = None,
        rna_emb: bool = False,
        zeroing_dapi: bool = False,
):
    """
    Parameters
    ----------
    epoch_index :
    tb_writer
    device
    rna2seg : model
    training_loader
    optimizer
    print_loss_every : print loss every x batch
    validation_loader : to compute validation loss (optinal)
    Returns
    -------
    """

    assert best_val_loss is not None
    assert cellbound_prob >= 0 and cellbound_prob <= 1, "cellbound_prob should be between 0 and 1"

    rna2seg.net.train()
    running_loss = 0.
    rna2seg.to(device)

    for i_train, dict_result in tqdm(enumerate(training_loader), file=sys.stdout,
                                     total=len(training_loader),
                                     desc="training"):  # file = sys.stdout, total=len(training_loader)):#

        loss = _run(dict_result=dict_result,
                    rna_emb=rna_emb,
                    device=device,
                    epoch_index=epoch_index,
                    no_cellbound_first_iter=no_cellbound_first_iter,
                    cellbound_index_to_use=cellbound_index_to_use,
                    cellbound_prob=cellbound_prob,
                    rna2seg=rna2seg,
                    optimizer=optimizer,
                    zeroing_dapi=zeroing_dapi)

        loss.backward()
        optimizer.step()

        running_loss += loss.item()
        if i_train % print_loss_every == 0:

            last_loss = running_loss / print_loss_every  # loss per batch

            tb_x = epoch_index * len(training_loader) + i_train + 1
            if tb_writer is not None:
                tb_writer.add_scalar('Loss/train', last_loss, tb_x)
            running_loss = 0.

            # validation loss

            if validation_loader is not None:
                rna2seg.eval()
                with torch.no_grad():
                    for i_test, dict_result_val in enumerate(validation_loader):
                        loss_val = _run(dict_result=dict_result_val,
                                        rna_emb=rna_emb,
                                        device=device,
                                        epoch_index=epoch_index,
                                        no_cellbound_first_iter=no_cellbound_first_iter,
                                        cellbound_index_to_use=cellbound_index_to_use,
                                        cellbound_prob=cellbound_prob,
                                        rna2seg=rna2seg,
                                        optimizer=optimizer,
                                        zeroing_dapi=zeroing_dapi,
                                        path_save_model=path_save_model,
                                        save_ouptut=len(training_loader) % max(i_train, 1) <= print_loss_every,
                                        i_test=i_test
                                        )

                        running_loss += loss_val.item()

                    # loss
                    last_loss_valid = running_loss / len(validation_loader)
                    print('  validation loss: {}'.format(last_loss_valid))
                    tb_x = epoch_index * len(training_loader) + i_train + 1
                    if tb_writer is not None:
                        tb_writer.add_scalar('Loss/valid', last_loss_valid, tb_x)

                    if best_val_loss is not None:
                        if last_loss_valid < best_val_loss:
                            best_val_loss = last_loss_valid
                            print(f"best_val_loss: {best_val_loss}")
                            torch.save(rna2seg.net.state_dict(),
                                       Path(path_save_model) / f"model_best_{epoch_index}_{i_train}.pt")
                            torch.save(rna2seg.net.state_dict(), Path(path_save_model) / f"best.pt")
                            if rna2seg.rna_embedding is not None:
                                torch.save(rna2seg.rna_embedding.state_dict(), Path(
                                    path_save_model) / f"model_best_rna_embedding{epoch_index}_{i_train}.pt")
                                torch.save(rna2seg.rna_embedding.state_dict(),
                                           Path(path_save_model) / f"best_rna_embedding.pt")

                    else:
                        best_val_loss = last_loss_valid
                        print(f"best_val_loss: {best_val_loss}")

                rna2seg.train()

                # save intermediate results

    return best_val_loss


def _concat_cellbound(dapi,
                      rna_img,
                      img_cellbound,
                      cellbound_prob=0.25,
                      cellbound_index_to_use=None,
                      zeroing_dapi=False,
                      ):
    assert dapi.ndim == 4
    assert img_cellbound.ndim == 4

    if not zeroing_dapi:
        for i in range(img_cellbound.shape[0]):
            if random.random() > cellbound_prob:
                img_cellbound[i] = 0
    else:
        for i in range(img_cellbound.shape[0]):
            if random.random() > cellbound_prob:
                img_cellbound[i] = 0
            if random.random() > cellbound_prob:
                dapi[i] = 0

    if cellbound_index_to_use is not None:
        img_cellbound = img_cellbound[:, cellbound_index_to_use[:img_cellbound.shape[1]]]

    img_concat = torch.cat((dapi, rna_img, img_cellbound), dim=1)
    return img_concat, dapi, rna_img, img_cellbound


def _run(dict_result,
         rna_emb,
         device,
         epoch_index,
         no_cellbound_first_iter,
         cellbound_index_to_use,
         cellbound_prob,
         rna2seg,
         optimizer,
         zeroing_dapi,
         path_save_model=None,
         save_ouptut=False,
         i_test=None):
    label, mask_gradient, img_cellbound, dapi, rna_img = (
        dict_result["mask_flow"], dict_result["mask_gradient"], dict_result["img_cellbound"],
        dict_result["dapi"], dict_result["rna_img"])
    if rna_emb:
        list_gene = dict_result['list_gene']
        array_coord = dict_result['array_coord']
        list_gene = list_gene.to(device)
        array_coord = array_coord.to(device)
    else:
        list_gene = None
        array_coord = None

    if epoch_index >= no_cellbound_first_iter:
        # sample by setinng to zero with a probability of cellbound_prob the img_cellbound
        imgs, dapi, rna_img, img_cellbound = _concat_cellbound(
            dapi=dapi,
            rna_img=rna_img,
            img_cellbound=img_cellbound,
            cellbound_prob=cellbound_prob,
            cellbound_index_to_use=cellbound_index_to_use,
            zeroing_dapi=zeroing_dapi,
        )
    else:  # no cellbound for the first iterations
        imgs, dapi, rna_img, img_cellbound = _concat_cellbound(dapi=dapi,
                                                               rna_img=rna_img,
                                                               img_cellbound=img_cellbound,
                                                               cellbound_prob=0,
                                                               cellbound_index_to_use=cellbound_index_to_use,
                                                               zeroing_dapi=zeroing_dapi)

    rna2seg.train()
    optimizer.zero_grad()
    label = label.to(device)
    dapi = dapi.to(device)
    img_cellbound = img_cellbound.to(device)
    rna_img = rna_img.to(device)

    dapi = dapi.to(torch.float32)
    rna_img = rna_img.to(torch.float32)
    img_cellbound = img_cellbound.to(torch.float32)

    output = rna2seg(
        list_gene=list_gene,
        array_coord=array_coord,
        dapi=dapi,
        img_cellbound=img_cellbound,
        rna_img=rna_img,
    )

    loss = partial_cellpose_loss(
        lbl=label[:, :3],
        y=output[:, :3],
        device=device,
        mask_gradient=mask_gradient
    )

    if path_save_model is not None and save_ouptut:  # (len(training_loader) % max(i_train, 1) <= print_loss_every):
        path_save_plot_epoch = Path(path_save_model) / f"epoch_{epoch_index}"
        path_save_plot_epoch.mkdir(exist_ok=True, parents=True)
        ##
        imgs = imgs.to("cpu")
        label = label.to("cpu")
        mask_gradient = mask_gradient.to("cpu")
        img_cellbound = img_cellbound.to("cpu")
        output = output.to("cpu")
        patch_index = dict_result["patch_index"]
        with open(path_save_plot_epoch / f"patch_index{i_test}.txt", "w") as f:
            f.write(str(patch_index.to("cpu").numpy()))

        for index in range(imgs.shape[0]):
            tifffile.imwrite(path_save_plot_epoch / f"imgs{patch_index[index]}.tif", imgs[index].numpy())
            tifffile.imwrite(path_save_plot_epoch / f"labels{patch_index[index]}.tif", label[index].numpy())
            tifffile.imwrite(path_save_plot_epoch / f"mask_gradient{patch_index[index]}.tif",
                             mask_gradient[index].numpy())
            tifffile.imwrite(path_save_plot_epoch / f"img_cellbounds{patch_index[index]}.tif",
                             img_cellbound[index].numpy())
            tifffile.imwrite(path_save_plot_epoch / f"output{patch_index[index]}.tif", output[index].numpy())

    return loss
