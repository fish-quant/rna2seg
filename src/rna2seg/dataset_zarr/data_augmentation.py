import albumentations as A
import cv2
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from cellpose.transforms import \
    random_rotate_and_resize as _random_rotate_and_resize


# Input image and Target image Transformations
def random_rotate_and_resize(input_image, mask_flow):
    shape = input_image.shape[-2]

    if mask_flow.shape[0] == 3:
        imgi, lbl, scale = _random_rotate_and_resize(
            [input_image], [mask_flow], xy=(shape, shape),
            rescale=[0.9, 0.9],  # why rescale is 0.9?
            random_per_image=False
        )
        imgi = imgi[0]
        lbl = lbl[0]

    elif mask_flow.shape[0] == 6:
        imgi, lbl, scale = _random_rotate_and_resize(
            [input_image, input_image],
            [mask_flow[:3], mask_flow[3:6]],
            xy=(shape, shape),
            rescale=[0.9, 0.9],
            random_per_image=False
        )
        imgi = imgi[0]
        lbl = np.concatenate((lbl[0], lbl[1]), axis=0)

    else:
        raise NotImplementedError("mask_flow shape not supported. n_channels should be 3 or 6.")
    assert imgi.shape == input_image.shape
    assert lbl.shape == mask_flow.shape

    return imgi, lbl, scale


# Cellbound stainings Transformations

class RandomChannelSwap(A.ImageOnlyTransform):
    def __init__(self, p=0.2, always_apply=False):
        super().__init__(p, always_apply)

    def apply(self, img, **params):
        if img.ndim == 3 and img.shape[2] >= 2:
            # "Image must have at least two channels to swap.", Skip if image has less than
            # 2 channels for dapi augmentation
            idx1, idx2 = np.random.choice(img.shape[2], 2, replace=False)
            img[..., [idx1, idx2]] = img[..., [idx2, idx1]]
        return img


class CoarseDropoutOneChannel(A.ImageOnlyTransform):
    def __init__(self, max_holes=8, max_height=8, max_width=8, min_holes=1, min_height=1, min_width=1, p=0.5):
        super().__init__(p)
        self.coarse_dropout = A.CoarseDropout(max_holes=max_holes, max_height=max_height, max_width=max_width,
                                              min_holes=min_holes, min_height=min_height, min_width=min_width, p=1.0)

    def apply(self, img, **params):
        if img.ndim == 3:
            channel = np.random.randint(0, img.shape[-1])
            img[..., channel] = self.coarse_dropout(image=img[..., channel].copy())['image']
        return img


cellbound_transform = A.Compose([
    A.RandomBrightnessContrast(brightness_limit=0.25, contrast_limit=0.25, p=0.2),
    A.GaussianBlur(blur_limit=(3, 7), p=0.2),
    A.CoarseDropout(max_holes=2, max_height=400, max_width=400, min_holes=1, min_height=1, min_width=1, p=0.2),
    CoarseDropoutOneChannel(max_holes=2, max_height=400, max_width=400, min_holes=1, min_height=100, min_width=100,
                            p=0.2),
    # A.HueSaturationValue(hue_shift_limit=0.01, sat_shift_limit=0.01, val_shift_limit=0.01, p=0.2),
    RandomChannelSwap(always_apply=False, p=0.2),
])


# RNA Embedding Transformations

def add_noise_to_gene_embeddings(gene2vect, noise_level=0.01, p=0.5):
    noisy_embeddings = {}
    for gene, emb in gene2vect.items():
        if np.random.rand() < p:
            noise = np.random.normal(loc=0.0, scale=noise_level, size=len(emb))
            noisy_embeddings[gene] = emb + noise
        else:
            noisy_embeddings[gene] = emb
    return noisy_embeddings


def shift_embeddings(gene2vect, min_scale_factor=0.9, max_scale_factor=1.1, min_shift=-1, max_shift=1, p=0.5):
    scaled_embeddings = {}
    for gene, emb in gene2vect.items():
        # independent shift of dimensions
        emb = np.array(emb)
        assert emb.ndim == 1, "Embedding must be 1D"

        for i in range(emb.shape[0]):
            if np.random.rand() < p:
                scale_factor = np.random.uniform(min_scale_factor, max_scale_factor)
                shift = np.random.uniform(min_shift, max_shift)  # Generate random shift
                scaled_embeddings[gene] = np.array(emb) * scale_factor + shift
                # emb[i] = emb[i] * np.random.uniform(min_scale_factor, max_scale_factor)
                # + np.random.uniform(min_shift, max_shift)
        scaled_embeddings[gene] = emb
    return scaled_embeddings


def augment_embeddings(gene2vect,
                       noise_level=0.01,
                       min_scale_factor=0.9,
                       max_scale_factor=1.1,
                       min_shift=-1,
                       max_shift=1,
                       noise_p=0.5,
                       shift_p=0.5):
    augmented = gene2vect.copy()
    augmented = add_noise_to_gene_embeddings(augmented, noise_level=noise_level, p=noise_p)
    augmented = shift_embeddings(
        augmented,
        min_scale_factor=min_scale_factor, max_scale_factor=max_scale_factor,
        min_shift=min_shift, max_shift=max_shift, p=shift_p)
    return augmented
