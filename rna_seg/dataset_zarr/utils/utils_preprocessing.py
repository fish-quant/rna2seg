


from tqdm import trange
import numpy as np
import os
import tifffile

def labels_to_flows_omnipose(labels, links=None, files=None, use_gpu=False, device=None,
                             omni=True, redo_flows=False, dim=2):
    """ Convert labels (list of masks or flows) to flows for training model.

    if files is not None, flows are saved to files to be reused

    Parameters
    --------------
    labels: list of ND-arrays
        labels[k] can be 2D or 3D, if [3 x Ly x Lx] then it is assumed that flows were precomputed.
        Otherwise labels[k][0] or labels[k] (if 2D) is used to create flows.
    links: list of label links
        These lists of label pairs define which labels are "linked",
        i.e. should be treated as part of the same object. This is how
        Omnipose handles internal/self-contact boundaries during training.
    files: list of strings
        list of file names for the base images that are appended with '_flows.tif' for saving.
    use_gpu: bool
        flag to use GPU for speedup. Note that Omnipose fixes some bugs that caused the Cellpose GPU
        implementation to have different behavior compared to the Cellpose CPU implementation.
    device: torch device
        what compute hardware to use to run the code (GPU VS CPU)
    omni: bool
        flag to generate Omnipose flows instead of Cellpose flows
    redo_flows: bool
        flag to overwrite existing flows. This is necessary when changing over from Cellpose to Omnipose,
        as the flows are very different.
    dim: int
        integer representing the intrinsic dimensionality of the data. This allows users to generate 3D flows
        for volumes. Some dependencies will need to be to be extended to allow for 4D, but the image and label
        loading is generalized to ND.

    Returns
    --------------
    flows: list of [4 x Ly x Lx] arrays
        flows[k][0] is labels[k], flows[k][1] is cell distance transform, flows[k][2:2+dim] are the
        (T)YX flow components, and flows[k][-1] is heat distribution / smooth distance

    """
    from omnipose.core import masks_to_flows

    nimg = len(labels)
    if links is None:
        links = [None]*nimg # just for entering below
    no_flow = labels[0].ndim != 3+dim # (6,Lt,Ly,Lx) for 3D, masks + dist + boundary + flow components, then image dimensions

    if no_flow or redo_flows:


        # compute flows; labels are fixed in masks_to_flows, so they need to be passed back
        with suppress_output():
            labels, dist, bd, heat, veci = map(list,zip(*[masks_to_flows(labels[n], links=links[n], use_gpu=use_gpu,
                                                                         device=device, omni=omni, dim=dim)
                                                          for n in trange(nimg)]))

        # concatenate labels, distance transform, vector flows, heat (boundary and mask are computed in augmentations)
        if omni:
            flows = [np.concatenate((labels[n][np.newaxis,:,:],
                                     dist[n][np.newaxis,:,:],
                                     veci[n],
                                     heat[n][np.newaxis,:,:]), axis=0).astype(np.float32)
                     for n in range(nimg)]
            # clean this up to swap heat and flows and simplify code? would have to rerun all flow generation
        else:
            flows = [np.concatenate((labels[n][np.newaxis,:,:],
                                     labels[n][np.newaxis,:,:]>0.5,
                                     veci[n]), axis=0).astype(np.float32)
                     for n in range(nimg)]
        if files is not None:
            for flow, file in zip(flows, files):
                file_name = os.path.splitext(file)[0]
                tifffile.imsave(file_name+'_flows.tif', flow)
    else:
        flows = [labels[n].astype(np.float32) for n in range(nimg)]

    return flows


### code to silent cellpose log



import sys
import os
from contextlib import contextmanager

@contextmanager
def suppress_output():
    # Save the current stdout and stderr
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    try:
        # Redirect stdout and stderr to /dev/null
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')
        yield
    finally:
        # Restore stdout and stderr
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

def noisy_function():
    print("This is a noisy function")