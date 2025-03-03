from __future__ import annotations


import torch
import torch.nn as nn
from einops import rearrange
from instanseg.utils.models.ChannelInvariantNet import (
    ChannelInvariantDecoderBlock, ChannelInvariantEncoderBlock)
from instanseg.utils.models.InstanSeg_UNet import conv_norm_act


class ChannelInvariantUNet(nn.Module):

    def __init__(self,
                 out_channels,
                 layers=[256, 128, 64, 32],
                 norm="BATCH",
                 act="ReLu",
                 aggregation="concat"):
        super().__init__()
        layers = layers[::-1]
        self.encoder = nn.ModuleList([ChannelInvariantEncoderBlock(1,
                                                                   layers[0],
                                                                   act=act,
                                                                   norm=norm,
                                                                   pool=False,
                                                                   aggregation=aggregation)] +

                                     [ChannelInvariantEncoderBlock(layers[i], layers[i + 1],
                                                                   act=act,
                                                                   norm=norm,
                                                                   pool=True,
                                                                   aggregation=aggregation)
                                      for i in range(len(layers) - 1)])
        layers = layers[::-1]
        self.decoder = nn.ModuleList([ChannelInvariantDecoderBlock(layers[i], layers[i + 1], layers[i + 1],
                                                                   act=act,
                                                                   norm=norm,
                                                                   final_decoder=(i == (len(layers) - 1)),
                                                                   aggregation=aggregation)
                                      for i in range(len(layers) - 1)])
        # finale conv
        final_norm = norm if (norm is not None) and norm.lower() != "instance" else None
        self.final_conv = conv_norm_act(layers[-1], out_channels, sz=1, norm=final_norm, act=None)

    def forward(self, x):

        b, c = x.shape[:2]
        x = rearrange(x, 'b c h w -> (b c) 1 h w')
        print(x.shape)

        # Encoder
        list_skip = []
        for n, layer in enumerate(self.encoder):
            x = layer(x, c=c, b=b)
            if n < len(self.encoder) - 1:
                list_skip.append(x)

        for n, layer in enumerate(self.decoder):
            x = layer(x, list_skip[-(n + 1)], c=c, b=b)

        x = self.final_conv(x)

        return x.float()


if __name__ == "__main__":
    ci_unet = ChannelInvariantUNet(out_channels=6,
                                   layers=[256, 128, 64, 32],
                                   norm="BATCH",
                                   act="ReLu",
                                   aggregation="concat")

    num_params = sum(p.numel() for p in ci_unet.parameters())
    print(f"Number of parameters in C_unet: {num_params}")

    print(ci_unet(torch.randn(1, 3, 600, 600)).shape)

    print(ci_unet(torch.randn(10, 3, 600, 600)).shape)
