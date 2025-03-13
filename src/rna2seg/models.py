
import dask

dask.config.set({'dataframe.query-planning': False})

import logging
import os
from pathlib import Path

import albumentations as A
import cv2
import geopandas as gpd
import numpy as np
import torch
import torch.nn as nn
from cellpose.dynamics import compute_masks
from cellpose.resnet_torch import CPnet
from instanseg.utils.models.ChannelInvariantNet import AdaptorNetWrapper
from sopa.segmentation import shapes
from torchvision.transforms.functional import gaussian_blur

log = logging.getLogger(__name__)


class CustomCPnet(CPnet):
    def __init__(self, *args, **kwargs):
        super(CustomCPnet, self).__init__(*args, **kwargs)

    def forward(self, imgs):
        t0 = self.downsample(imgs)
        style = self.make_style(t0[-1])
        t0 = self.upsample(style, t0, mkldnn=False)
        out = self.output(t0)
        return out


try:
    from vmunet.vmunet import VMUNet


    class CustomVMUnet(VMUNet):
        def __init__(self, *args, **kwargs):
            super(CustomVMUnet, self).__init__(*args, **kwargs)

        def load_model(self, filename):
            model_dict = self.vmunet.state_dict()
            pretrained_dict = torch.load(filename)
            pretrained_dict = {k.split("vmunet.")[1]: v for k, v in pretrained_dict.items()}
            new_dict = {
                k: v for k, v in pretrained_dict.items()
                if k in model_dict.keys()
            }
            model_dict.update(new_dict)
            print('Total model_dict: {}, Total pretrained_dict: {}, update: {}'.format(len(model_dict),
                                                                                       len(pretrained_dict),
                                                                                       len(new_dict)))
            self.vmunet.load_state_dict(model_dict)

            not_loaded_keys = [k for k in pretrained_dict.keys() if k not in new_dict.keys()]
            print('Not loaded keys:', not_loaded_keys)
            print("Encoder loaded finished!")

except Exception as e:
    print(e)
    print("VMUnet not loaded")


class RNAEmbedding(nn.Module):
    def __init__(
            self,
            gene2index,
            embedding_dim,
            special_index=0,
            device=None,
            gaussian_kernel_size=0,
            sigma=None,
            radius_rna=None,
    ):

        super(RNAEmbedding, self).__init__()
        self.gene2index = gene2index
        self.n_genes = len(self.gene2index) + 1  # add one for the special index
        self.embedding_dim = embedding_dim
        self.embedding = nn.Embedding(self.n_genes, embedding_dim)
        self.special_index = special_index
        # Initialize the special index with zeros and make it non-trainable
        self.embedding.weight.data[self.special_index] = torch.zeros(self.embedding_dim)
        self.gaussian_kernel_size = gaussian_kernel_size
        self.device = device
        self.sigma = sigma
        self.radius_rna = radius_rna
        self.to(device)

    def forward(self, shape, list_gene, array_coord
                ):

        # convert gene to index

        list_y = array_coord[:, :, 0]
        list_x = array_coord[:, :, 1]

        list_y = list_y.long()
        list_x = list_x.long()

        # create tensor of zero of  shape (batch_size, embedding_dim, img.shape[-2], img.shape[-1])
        rna_imgs = torch.zeros(shape[0], shape[-2], shape[-1], self.embedding_dim, device=self.device)
        emb = self.embedding(list_gene)
        assert torch.sum(emb[list_gene == self.special_index]) == 0
        emb[list_gene == self.special_index] = torch.zeros(self.embedding_dim,
                                                           device=self.device)

        batch_size = shape[0]
        # Vectorized approach:
        if self.radius_rna is None:
            # Use scatter operation for all batches at once
            for b in range(batch_size):
                rna_imgs[b].index_put_((list_y[b], list_x[b]), emb[b])
        else:
            # For the radius case, precompute offsets
            offsets_y, offsets_x = torch.meshgrid(
                torch.arange(-self.radius_rna, self.radius_rna + 1, device=self.device),
                torch.arange(-self.radius_rna, self.radius_rna + 1, device=self.device)
            )
            offsets = torch.stack([offsets_y.flatten(), offsets_x.flatten()], dim=1)
            
            for b in range(batch_size):
                for offset_y, offset_x in offsets:
                    y_coords = torch.clamp(list_y[b] + offset_y, 0, shape[-2] - 1)
                    x_coords = torch.clamp(list_x[b] + offset_x, 0, shape[-1] - 1)
                    rna_imgs[b].index_put_((y_coords, x_coords), emb[b])

        rna_imgs = rna_imgs.permute(0, 3, 1, 2)

        # add max filter with scipy
        if self.gaussian_kernel_size > 0:
            # todo : add a parameter to choose the sigma or kernel_size
            # todo  :  choose only one parameter to be set
            rna_imgs = gaussian_blur(rna_imgs,
                                     kernel_size=self.gaussian_kernel_size,
                                     sigma=self.sigma)

        return rna_imgs


class RNA2seg(nn.Module):
    """
    RNA2seg: A deep learning-based method for cell segmentation using 
    spatial transcriptomics data and membrane and nuclei stainings.
    """

    def __init__(
            self,
            device,
            pretrained_model: Path | str | None = "default_pretrained",
            net: str = "unet",
            n_inv_chan: int = 3,
            nb_rna_channel: int = 1,
            nout: int = 3,
            nbase=[32, 64, 128, 256],
            sz: int = 3,
            # diameter: int = 30,
            flow_threshold: float = 0.9,
            min_cell_size: float = 200,
            cellbound_flow_threshold: float = 0.4,
            gene2index=None,
    ):

        """
        Initialize RNA2Seg.

        :param device: The computing device (e.g., "cpu" or "cuda").
        :type device: str
        :param pretrained_model: Path to a pretrained model.
        If "default_pretrained", a trained rna2seg model is download from huggingface. If None weight are randomly
        initialized,  Defaults is "default_pretrained".
        :type pretrained_model: Path | str | None
        :param net: Backbone network architecture. Can be "unet" or "vmunet". Defaults to "unet".
        :type net: str
        :param n_inv_chan: Number of channels in staining input image. The stainings are combined and encoded into an
        image of n_inv_chan channels using a Channel-Net. Defaults to 3.
        :type n_inv_chan: int
        :param nb_rna_channel: Number of RNA channels used as input. Defaults to 1.
        :type nb_rna_channel: int
        :param nout: Number of output channels. Following the CellPose method, the network outputs cell probabilities on
         one channel and the 2-channel flow. Defaults to 3.
        :type nout: int
        :param nbase: List defining the number of channels at each layer of the network.
        :type nbase: list[int]
        :param sz: Kernel size for convolutions. Defaults to 3.
        :type sz: int
        :param flow_threshold: Threshold for flow consistency during segmentation. Defaults to 0.9.
        :type flow_threshold: float
        :param min_cell_size: Minimum cell size (in pixels) to retain. Defaults to 200.
        :type min_cell_size: float
        :param cellbound_flow_threshold: Threshold for cell boundary detection. Defaults to 0.4.
        :type cellbound_flow_threshold: float
        :param gene2index: Mapping from gene names to indices, if applicable. Defaults to None.
        :type gene2index: dict or None
        :raise ValueError: If an invalid model or configuration is provided.
        """
        super().__init__()
        assert nb_rna_channel is not None, "rna_channel must be provided"
        self.nb_rna_channel = nb_rna_channel

        self.device = device
        self.nout = nout
        self.sz = sz
        self.flow_threshold = flow_threshold
        self.min_cell_size = min_cell_size
        self.pretrained_model = pretrained_model

        self.cellbound_flow_threshold = cellbound_flow_threshold

        self.net_archi = net
        self.n_inv_chan = n_inv_chan

        if gene2index is not None:
            assert 0 not in gene2index.values(), "gene2index should not contain 0 as value"
            log.info("initilization RNA embedding layer")
            print("initilization RNA embedding layer")
            self.embedding_dim = 3
            self.rna_embedding = RNAEmbedding(
                gene2index=gene2index,
                embedding_dim=self.embedding_dim,
                device=self.device,
                gaussian_kernel_size=3,
                sigma=0.5,
                radius_rna=1,
            )
            self.rna_embedding = self.rna_embedding.to(device)
            self.nb_rna_channel = self.embedding_dim
        else:
            self.rna_embedding = None
            assert nb_rna_channel is not None, "rna_channel must be provided"
            self.nb_rna_channel = nb_rna_channel

        self._set_net(nbase)
        self.net = self.net.to(device)

    def forward(
            self,
            input_dict=None,
            list_gene=None,
            array_coord=None,
            dapi=None,
            img_cellbound=None,
            rna_img=None,
    ):  # to modify :  clean how RNA and staining are put in the image

        """
        Forward pass for the RNA2seg model.

        The `forward` method supports two modes of input:
            1. You can pass a dictionary (`input_dict`) containing all the relevant inputs ('list_gene', 'array_coord',
            'dapi', 'img_cellbound', and 'rna_img'). If provided, the values in the dictionary will override
            the individual function arguments.
            2. Alternatively, you can pass each argument independently.

        The DAPI and cell boundary images are used for encoding, and the RNA image (either encoded or pre-encoded) is
         combined with the other inputs. The concatenated data is then passed through the model to generate the output.
    
        
        :param input_dict: Optional dictionary containing the following keys: 'list_gene', 'array_coord', 'dapi',
        'img_cellbound', and 'rna_img'. If provided, the values will override the function arguments.
        :type input_dict: dict | None
        :param list_gene: Tensor representing the list of genes for RNA encoding. Cannot be provided simultaneously with
        `rna_img`. Defaults to None.
        :type list_gene: torch.Tensor | None
        :param array_coord: Tensor containing the coordinates for RNA encoding. Required if `list_gene` is provided.
        :type array_coord: torch.Tensor | None
        :param dapi: Tensor representing the DAPI staining image, used as input for encoding.
        :type dapi: torch.Tensor
        :param img_cellbound: Tensor representing the cell boundary image, used as input for encoding.
        :type img_cellbound: torch.Tensor
        :param rna_img: Tensor representing the pre-encoded RNA image. Cannot be provided simultaneously
         with `list_gene`. Defaults to None.
        :type rna_img: torch.Tensor | None

        :returns: Tensor representing the model's output after processing the input data.
        :rtype: torch.Tensor

        :raises AssertionError: If neither `list_gene` nor `rna_img` is provided, or if both are provided simultaneously.
        :raises AssertionError: If `array_coord` is not provided when `list_gene` is used.
        """

        if input_dict is not None:
            list_gene = input_dict.get('list_gene', list_gene)
            array_coord = input_dict.get('array_coord', array_coord)
            dapi = input_dict.get('dapi', dapi)
            img_cellbound = input_dict.get('img_cellbound', img_cellbound)
            rna_img = input_dict.get('rna_img', rna_img)

        dapi = dapi.to(self.device)
        img_cellbound = img_cellbound.to(self.device)
        if rna_img is not None:
            rna_img = rna_img.to(self.device)
        assert list_gene is not None or rna_img is not None, "list_gene or rna_img must be provided"
        assert list_gene is None or rna_img is None, \
            "list_gene to encode rna and rna_img using pre-encoded RNA cannot be provided at the same time"

        if list_gene is not None:  # encode RNA

            assert array_coord is not None, "array_coord must be provided if list_gene is provided to encode the RNA"

            list_gene = list_gene.to(self.device)
            array_coord = array_coord.to(self.device)
            rna_img = self.rna_embedding(shape=dapi.shape,
                                         list_gene=list_gene,
                                         array_coord=array_coord, )

        x = self.encode(dapi=dapi,
                        img_cellbound=img_cellbound)

        x = torch.cat((x, rna_img), dim=1)
        out = self.net.model(x)

        return out

    def encode(
            self,
            imgs=None,
            dapi=None,
            img_cellbound=None
    ):  # to modify :  clean how RNA and staining are put in the image

        """
        Encodes the input images (DAPI and cell boundary) for the model.

        The method supports two modes of input:
            1. You can pass a tensor `imgs` containing the full image data.
            2. Alternatively, you can pass DAPI and cell boundary images separately (`dapi` and `img_cellbound`).
             These two images are concatenated along the channel dimension and then passed through the network.
        
        :param imgs: Tensor containing the full set of image channels, including RNA channels.
        This is kept for compatibility with the old version. Defaults to None.
        :type imgs: torch.Tensor | None
        :param dapi: Tensor representing the DAPI staining image, used for encoding. Must be provided if `img_cellbound`
         is provided.
        :type dapi: torch.Tensor | None
        :param img_cellbound: Tensor representing the cell boundary image, used for encoding. Must be provided if `dapi`
         is provided.
        :type img_cellbound: torch.Tensor | None

        :returns: Tensor representing the encoded image data, after processing through the network's `AdaptorNet`.
        :rtype: torch.Tensor

        :raises AssertionError: If `imgs` is provided while `dapi` and `img_cellbound` are also provided, or if neither
         is provided when one of them is required.
        """

        if imgs is not None:  # old version to delete, only kept for compatibility
            assert dapi is None
            assert img_cellbound is None
            list_channel_to_merge = [0] + list(range(1 + self.nb_rna_channel, imgs.shape[1]))
            out = self.net.model.AdaptorNet(imgs[:, list_channel_to_merge, :, :])

        else:  # new version
            assert imgs is None
            if dapi is not None:
                assert img_cellbound is not None
                merge_staining = torch.cat([dapi, img_cellbound], dim=1)
                out = self.net.model.AdaptorNet(merge_staining)

        return out

    def run(
            self,
            path_temp_save,
            input_dict=None,
            list_gene=None,
            array_coord=None,
            dapi=None,
            img_cellbound=None,
            rna_img=None,
            bounds=None,
            min_area=0,
    ):

        """
        Evaluates the model on a batch of images or a single image, and optionally on staining images.

        :param path_temp_save: The directory where the results will be saved.
        :type path_temp_save: str | Path
        :param input_dict: A dictionary containing the inputs for the model. It can include 'list_gene', 
                    'array_coord', 'dapi', 'img_cellbound', 'rna_img', and 'bounds'.
        :type input_dict: dict | None
        :param list_gene: List of gene expressions to use for encoding RNA. Defaults to None.
        :type list_gene: torch.Tensor | None
        :param array_coord: Coordinates array for the genes, required if `list_gene` is provided.
        :type array_coord: torch.Tensor | None
        :param dapi: DAPI stained image used for encoding.
        :type dapi: torch.Tensor
        :param img_cellbound: Image of cell boundaries used for encoding.
        :type img_cellbound: torch.Tensor
        :param rna_img: RNA image, either encoded or pre-encoded.
        :type rna_img: torch.Tensor | None
        :param bounds: Bounds for the image, used for transformations. Defaults to None.
        :type bounds: list | None
        :param min_area: The minimum area to consider for detected cells. Defaults to 0 (no filtering).
        :type min_area: int

        :return: A tuple containing the flow, cell probability, predicted masks, and cells (as a GeoDataFrame).
        :rtype: tuple (torch.Tensor, torch.Tensor, np.array, GeoDataFrame)
        """
        Path(path_temp_save).mkdir(parents=True, exist_ok=True)
        if input_dict is not None:
            list_gene = input_dict.get('list_gene', list_gene)
            array_coord = input_dict.get('array_coord', array_coord)
            dapi = input_dict.get('dapi', dapi)
            img_cellbound = input_dict.get('img_cellbound', img_cellbound)
            rna_img = input_dict.get('rna_img', rna_img)
            bounds = input_dict.get('bounds', bounds)

        self.net.eval()
        if "rna_embedding" in self.__dict__ and self.rna_embedding is not None:
            self.rna_embedding.eval()

        # Check if the input is a batch or a single image
        if dapi.dim() == 3:
            dapi = dapi.unsqueeze(0)
            img_cellbound = img_cellbound.unsqueeze(0)
            if list_gene is not None:
                list_gene = list_gene.unsqueeze(0)
                array_coord = array_coord.unsqueeze(0)
            else:
                rna_img = rna_img.unsqueeze(0)


        # Forward
        res = self.forward(list_gene=list_gene,
                           array_coord=array_coord,
                           dapi=dapi,
                           img_cellbound=img_cellbound,
                           rna_img=rna_img)

        # Compute masks
        flow_list = []
        cellprob_list = []
        masks_pred_list = []
        for k in range(len(res)):
            flow, cellprob = res[k][:2], res[k][2]
            masks_pred = compute_masks(
                dP=flow.to("cpu").long().detach().numpy(),
                cellprob=cellprob.to("cpu").long().detach().numpy(),
                min_size=self.min_cell_size,
                flow_threshold=self.flow_threshold,
                cellprob_threshold=-0.5,
                niter=200, interp=True, do_3D=False,
                device=flow.device,
            )
            flow_list.append(flow)
            cellprob_list.append(cellprob)
            if type(masks_pred) == tuple:
                masks_pred = masks_pred[0]
            masks_pred_list.append(masks_pred)

        flow = torch.stack(flow_list, dim=0)
        cellprob = torch.stack(cellprob_list, dim=0)
        masks_pred = np.stack(masks_pred_list, axis=0)

        # Transform to original size
        assert input_dict['bounds'][0] - input_dict['bounds'][2] == input_dict['bounds'][1] - input_dict['bounds'][
            3], "image must be square"
        original_image_shape = input_dict['bounds'][2] - input_dict['bounds'][0]
        transforms_img1 = A.Compose([
            A.Resize(
                width=original_image_shape,
                height=original_image_shape,
                interpolation=cv2.INTER_NEAREST
            ),
        ])
        masks_pred = transforms_img1(image=masks_pred[0])["image"]

        # Create cells shape
        cells = shapes.geometrize(masks_pred, tolerance=None, smooth_radius_ratio=0.1)
        print(f'{len(cells)} cells detected')
        cells = cells[cells.area >= self.min_area] if min_area > 0 else cells
        cells = gpd.GeoDataFrame(geometry=cells.geometry)
        cells.geometry = cells.translate(*input_dict['bounds'][:2])

        # Save the cells as parquet
        if path_temp_save is not None and input_dict is not None:
            cells.to_parquet(path_temp_save / f"{input_dict['idx']}.parquet")

        return flow, cellprob, masks_pred, cells

    def save_model(self, filename):
        """
        Save the model to a file.

        :param filename: The path to the file where the model will be saved.
        :type filename: str
        """
        torch.save(self.net.state_dict(), filename)

    def load_model(self, filename, device=None):

        """
        Load the model from a file.

        :param filename: The path to the file where the model is saved.
        :type filename: str
        :param device: The device to load the model on. If None, the model is loaded on the CPU. 
                    Defaults to None.
        :type device: torch.device | None
        """
        model_dict = torch.load(filename, map_location='cpu')
        self.net.load_state_dict(model_dict, strict=True)
        self.net.to(device)

    def _set_net(self, nbase, ):

        nchan = self.n_inv_chan + self.nb_rna_channel
        if self.net_archi == "unet":
            print("initiaisation of CPnet")
            nbase = [nchan, *nbase]
            self.net = CustomCPnet(nbase=nbase, nout=self.nout, sz=self.sz, mkldnn=False) #, diam_mean=self.diameter)
        elif self.net_archi == "vmunet":
            print("initiaisation of VMUNet")
            self.net = CustomVMUnet(num_classes=self.nout, input_channels=nchan)
        else:
            raise ValueError(f"Model not implemented: {self.net}")

        print("Initiaisation of ChannelInvariantNet")
        model = self.net
        self.net = AdaptorNetWrapper(model, out_channels=self.n_inv_chan)

        self.net = self.net.to(self.device)
        if self.pretrained_model == "default_pretrained":
            from huggingface_hub import snapshot_download
            path = snapshot_download(repo_id="aliceblondel/RNA2seg")
            print(f"Downloading pretrained model from huggingface {path}")
            self.pretrained_model = Path(path) / "rna2seg.pt"

        if self.pretrained_model == "brain_hamster":
            from huggingface_hub import snapshot_download
            path = snapshot_download(repo_id="faraone75/rna2seg_hamster_brain")
            print(f"Downloading pretrained model from huggingface {path}")
            self.pretrained_model = Path(path) / "brain_hamster.pt"

        if self.pretrained_model is not None:
            assert os.path.exists(self.pretrained_model), f"Pretrained model not found at : {self.pretrained_model}"
            print(f"Loading weights from {self.pretrained_model}")
            self.load_model(self.pretrained_model, device=self.device)
