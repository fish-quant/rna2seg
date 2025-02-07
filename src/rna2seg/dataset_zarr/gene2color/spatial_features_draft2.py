




#%%
import time
import numpy as np
import os
import pandas
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import tifffile
from sklearn.neighbors import NearestNeighbors
import networkx as nx
import scipy.sparse as sp
import scipy
from tqdm import tqdm
from pathlib import Path


## import ndi
from scipy import ndimage as ndi
from sklearn.neighbors import NearestNeighbors


def compute_rna2nuclei_distance(nuclei_all_z,
                                df_spots,
                                plot_distribution = False):
    """

    Args:
        nuclei_all_z:  np.array of shape (z, y, x) containing the mask of the nuclei for each z
        df_spots: dataframe containing the spots with the columns x, y, z, gene
        plot_distribution:

    Returns:
        compute the "relative" mean and std distance of eahc RNA species to the nuclei boundary
        relative mean here negative if inside, positive if outside
        distance2nuclei_boundary_mean is a pandas series with the gene as index and the mean distance as value
        distance2nuclei_boundary_std is a pandas series with the gene as index and the std distance as value

    """

    try_contour3D = np.zeros(nuclei_all_z.shape)
    for z in tqdm(list(range(nuclei_all_z.shape[0]))):
        try_contour3D[z] = ndi.maximum_filter(nuclei_all_z[z], size=2) - ndi.minimum_filter(
            nuclei_all_z[z], size=2)
    try_contour3D_bool = try_contour3D != 0

    distance_all_z = np.zeros(nuclei_all_z.shape)
    for z in tqdm(list(range(nuclei_all_z.shape[0]))):
        inverted_mask = np.ones(try_contour3D_bool[z].shape, dtype = np.uint8) - (try_contour3D_bool[z] != 0).astype(np.uint8)
        distance = ndi.distance_transform_edt(inverted_mask,
                                              sampling=[0.108, 0.108])
        distance_all_z[z] = distance
    np.save("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/distance_contour3D.npy", distance_all_z)

    f_sign = lambda x: -1 if x <= 0 else 1
    list_distance2nuclei_boundary = []
    list_is_insed_nuclei = []
    for index, spots in tqdm(df_spots.iterrows(), total = len(df_spots)):
        x = spots.x
        y = spots.y
        z = spots.z
        distance2nuclei_boundary = f_sign(nuclei_all_z[z, y, x]) * distance_all_z[z, y, x]
        list_distance2nuclei_boundary.append(distance2nuclei_boundary)
        list_is_insed_nuclei.append(f_sign(nuclei_all_z[z, y, x]) == -1)
    df_spots["distance2nuclei_boundary"] = list_distance2nuclei_boundary
    df_spots["is_inside_nuclei"] = list_is_insed_nuclei
    ### compute the mean and std of the distance2nuclei_boundary for each gene with a group by gene
    distance2nuclei_boundary_mean = df_spots.groupby("gene")["distance2nuclei_boundary"].mean()
    distance2nuclei_boundary_std = df_spots.groupby("gene")["distance2nuclei_boundary"].std()

    if plot_distribution:
        for gene in distance2nuclei_boundary_mean.index:
            sns.histplot(df_spots[df_spots.gene == gene].distance2nuclei_boundary, bins=100, kde=True)
            plt.title(gene)
            plt.show()
    return distance2nuclei_boundary_mean, distance2nuclei_boundary_std, df_spots





#%%





if __name__ == "__main__":


    nuclei_all_z = tifffile.imread("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/mask_nuclei_before_erosion/nuclei.tif")
    #df_spots = pd.read_csv("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/csv/cellpose_dapi.csv")
    df_spots = pd.read_csv("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/csv/prior_stitched.csv")
    gene2coexpression = np.load("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/gene_color/gene2color_npc59.npy", allow_pickle=True).item()


    distance2nuclei_boundary_mean, distance2nuclei_boundary_std, df_spots = compute_rna2nuclei_distance(nuclei_all_z,
                                                                                              df_spots,
                                                                                              plot_distribution=False)


    dict_distance2nuclei_boundary_mean = distance2nuclei_boundary_mean.to_dict()
    ## sort gene by distance2nuclei_boundary_mean
    sorted_gene = sorted(dict_distance2nuclei_boundary_mean, key=dict_distance2nuclei_boundary_mean.get)


    dict_distance2nuclei_boundary_std = distance2nuclei_boundary_std.to_dict()
    ## sort gene by distance2nuclei_boundary_std
    sorted_gene = sorted(dict_distance2nuclei_boundary_std, key=dict_distance2nuclei_boundary_std.get)

    dict_distance2nuclei_boundary_mean['Neat1']
    ### gene_list = ['Apob', 'Neat1', "Slc5a1"]


    gene_list_to_plot = ['Apob', 'Neat1', "Slc5a1"]
    dict_color = {"Apob": "orange", "Neat1": "blue", "Slc5a1": "green"}


    df_spots.to_csv("/home/tom/Bureau/FIGURE/sympAI_pasteur/mouse_ileum/rna_loc.csv")
    plt.rcParams['xtick.labelsize'] = 45
    plt.rcParams['ytick.labelsize'] = 45
    plt.rc_context({"axes.labelsize": 25, })

    fig , ax = plt.subplots(1, 1, figsize = (13, 10))
    for gene in  ['Apob', "Neat1"]:
        list_mean  = list(df_spots[df_spots.gene == gene].distance2nuclei_boundary)
        print(len(list_mean))
        ## sample 10000 points
        list_mean = np.random.choice(list_mean, 10000)
        #sns.kdeplot(list_mean, color=dict_color[gene], label=gene, ax=ax,
         #linewidth=5)
        sns.histplot(list_mean,
                     bins=150, kde=True, color = dict_color[gene],
                        label = gene,ax=ax, linewidth=5)
    ax.set_xlim(-4, 4)
    # #ax.set_yscale('log')
    #ax.legend()
    plt.show()

    plot_title = "distance2nuclei_boundary"
    folder_save = '/home/tom/Bureau/FIGURE/sympAI_pasteur/mouse_ileum/rna_loc/'
    image_format = 'svg' # e.g .png, .svg, etc.
    file_name = Path(folder_save) / f'{plot_title}.svg'
    fig.savefig(file_name, format=image_format, dpi=200)
    image_format = 'png'  # e.g .png, .svg, etc.
    file_name = Path(folder_save) / f'{plot_title}.png'
    fig.savefig(file_name, format=image_format, dpi=200)
    plt.show()





    ### plot the distribution of the distance2nuclei_boundary for each gene


    df_spots_label = compute_coexpression_heterogeneity_per_spots(
        df_spots_label = df_spots,
        gene2coexpression = gene2coexpression,
        dict_scale={"x": 0.108, "y": 0.108, "z": 1},
        z_column_name="z",
        y_column_name="y",
        x_column_name="x",
        radius=3,
        n_neighbors=15,
        remove_self_node=True,
        list_z=list(range(2, 6)))

    type(distance2nuclei_boundary_std)

    ### plot the distribution of the distance2nuclei_boundary for each gene
    for gene in distance2nuclei_boundary_mean.index:



    #### plot vizualization of the distance2nuclei_boundary for each gene

    ### generate a color vector for each gene
    assert list(distance2nuclei_boundary_std.index)  == list(distance2nuclei_boundary_mean.index)
    norm_distance2nuclei_boundary_mean = np.array(distance2nuclei_boundary_mean) - np.array(distance2nuclei_boundary_mean).min()

    norm_distance2nuclei_boundary_mean /= np.array(norm_distance2nuclei_boundary_mean).max()

    cmap = plt.get_cmap('cividis')

    cmap = plt.get_cmap('viridis')
    dict_color = {gene: np.array(list(cmap(value)[:3])) for gene, value in zip(distance2nuclei_boundary_mean.index, norm_distance2nuclei_boundary_mean)}
    np.save("/home/tom/Bureau/FIGURE/sympAI_pasteur/mouse_ileum/dict_color_rgb_dist.npy", dict_color)

    dict_color = {}
    for gene in distance2nuclei_boundary_mean.index:
        gene_index = list(list(distance2nuclei_boundary_std.index)).index(gene)
        dict_color[gene] = [1 - norm_distance2nuclei_boundary_mean[gene_index],
                            norm_distance2nuclei_boundary_mean[gene_index],
                            1 - norm_distance2nuclei_boundary_mean[gene_index]]


    from matplotlib import colors

    # Generate a color vector for each gene
    assert list(distance2nuclei_boundary_std.index) == list(distance2nuclei_boundary_mean.index)
    norm_distance2nuclei_boundary_mean = np.array(distance2nuclei_boundary_mean) - np.array(distance2nuclei_boundary_mean).min()
    norm_distance2nuclei_boundary_mean /= np.array(norm_distance2nuclei_boundary_mean).max()


    gene_list = ['Apob', 'Neat1', "Slc5a1"]# + list(df_spots.gene.unique())[:50]
    import napari
    dapi = tifffile.imread(
        "/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/images/nuc_dapi_z4.tif")
    df2 = pd.read_csv("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/rna_z_stack/z_3.csv")
    df3 = pd.read_csv("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/rna_z_stack/z_4.csv")
    df4 = pd.read_csv("/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/rna_z_stack/z_5.csv")
    df = pd.concat([df2, df3, df4])

    #df.global_y = df.global_y / 0.108
    #df.global_x = df.global_x / 0.108


    viewer = napari.Viewer()
    gene_list = ['Apob', 'Neat1', "Slc5a1"] + sorted_gene[:25] + sorted_gene[-25:]
    viewer.add_image(dapi, name='dapi')
    for gene in tqdm(gene_list[:]):
        rna_gene = df[df["gene"] == gene]
        ## random color
        # color = gene_color_dico[gene]
        #color = "#" + "%06x" % np.random.randint(0, 0xFFFFFF)
        color = dict_color[gene]
        viewer.add_points(rna_gene[["y", "x"]], name=gene, face_color=color, size=9)

    viewer.add_image(dapi, name='dapi')

#/home/tom/Bureau/FIGURE/sympAI_pasteur/mouse_ileum/


    viewer = napari.Viewer()

    viewer.add_image(image, name='dapi')
    viewer.add_image(cyto1, name='cyto1')
    viewer.add_image(cyto2, name='cyto2')
    viewer.add_image(cyto3, name='cyto3')
    viewer.add_image(rna, name='rna')





    ### plot coexpression heterogeneity


    dapi = tifffile.imread(
        "/media/tom/Transcend1/mouse_ileum_rnaseg/crop1/images/nuc_dapi_z4.tif")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.imshow(dapi,
              cmap="gray")
    plt.show()

    df_spots_label_z345 = df_spots_label[np.isin(df_spots_label.z,
                                                 [4])].copy()

    #norm min max
    df_spots_label_z345.coexpression_heterogeneity = ((df_spots_label_z345.coexpression_heterogeneity - df_spots_label_z345.coexpression_heterogeneity.min()) /
                                                      (df_spots_label_z345.coexpression_heterogeneity.max() - df_spots_label_z345.coexpression_heterogeneity.min()))

    y_min = 3800
    y_max = 4500
    x_min = 1800
    x_max = 4000
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.imshow(dapi[y_min:y_max, x_min:x_max],
              cmap="gray")
    spot_size = 4
    for index, spots in tqdm(df_spots_label_z345.iterrows(), total=len(df_spots_label_z345)):
        if y_min < spots.y < y_max and x_min < spots.x < x_max:
            x = spots.x
            y = spots.y
            z = spots.z
            color = [1 - spots.coexpression_heterogeneity, spots.coexpression_heterogeneity,
                     1 - spots.coexpression_heterogeneity]
            ax.scatter( x-x_min, y-y_min, color=color, s=spot_size)
    plt.show()