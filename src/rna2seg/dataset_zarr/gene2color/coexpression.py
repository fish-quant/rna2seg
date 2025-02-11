

#%%
import numpy as np
import pandas
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import scipy.sparse as sp
import scipy
from tqdm import tqdm
from pathlib import Path
######### use the coexpression matrix to embedd the gene in 2D where similar gene are close to each other
from rna2seg._constant import RNA2segFiles


def count_matrix_in_situ_from_knn(list_gene,
                                  df_spots_label,
                                  n_neighbors=5,
                                  radius=None,
                                  remove_self_node = False,
                                  dict_scale = {"x": 0.108, "y": 0.108, "z": 1},
                                  sampling=True,
                                  sampling_size=10000,
                                  z_column_name = "z",
                                  y_column_name = "y",
                                  x_column_name = "x"
                                  ):


    """
    Compute the colon-expression score matrix for the RNA spatial distribution

    :param df_spots_label:  dataframe with the columns x,y,z,gene. the coordinates are rescaled in µm by dict_scale attribute of the dataset object
    :type df_spots_label: pd.DataFrame
    :param n_neighbors: maximum number of neighbors default is 40
    :type n_neighbors: int
    :param radius: maximum radius of neighbors. It should be set proportionnaly to expected cell size, default is radius =  mean_cell_diameter / 2
    :return: count_matrix of shape (N_rna,  n_genes) where n_genes is the number of unique genes in df_spots_label
    each row is an 'RNA expression vector' summarizing local expression neighborhood of a molecule
    :rtype: np.array
    """

    gene_index_dico = {}
    for gene_id in range(len(list_gene)):
        gene_index_dico[list_gene[gene_id]] = gene_id #todo gene_index_dico to add in self
    ## this part should be factorize with create_directed_nn_graph
    try:
        df_spots_label = df_spots_label.reset_index()
    except Exception as e:
        print(e)
    if z_column_name in df_spots_label.columns:
        list_coordo_order_no_scaling = np.array([df_spots_label[z_column_name], df_spots_label[y_column_name], df_spots_label[x_column_name]]).T
        list_coordo_order = list_coordo_order_no_scaling * np.array(
            [dict_scale['z'], dict_scale['y'],dict_scale["x"]])
    else:
        list_coordo_order_no_scaling = np.array([df_spots_label[y_column_name], df_spots_label[x_column_name]]).T
        list_coordo_order = list_coordo_order_no_scaling * np.array([dict_scale['y'], dict_scale['x']])

    dict_list_features = {}
    assert 'gene' in df_spots_label.columns
    for feature in df_spots_label.columns:
        dict_list_features[feature] = list(df_spots_label[feature])
    array_gene_indexed = np.array([dict_list_features['gene'][i] for i in range(len(df_spots_label))])


    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(list_coordo_order)
    ad = nbrs.kneighbors_graph(list_coordo_order)  ## can be optimize here
    distance = nbrs.kneighbors_graph(list_coordo_order, mode='distance')
    ad[distance > radius] = 0
    ad.eliminate_zeros()


    rows, cols, BOL = sp.find(ad == 1)

    unique_rows = np.unique(rows)
    if sampling:
        if len(unique_rows) > sampling_size:
            unique_rows = np.random.choice(unique_rows, sampling_size, replace=False)

    list_expression_vec = []
    for row in tqdm(unique_rows):
        col_index = np.nonzero(ad[row])[1]
        if remove_self_node:
            col_index = col_index[col_index != row]
        vectors_gene = list(array_gene_indexed[col_index])
        vector_distance = np.array([distance[row, col] for col in col_index])
        expression_vector = np.zeros(len(list_gene))
        for str_gene_index in range(len(vectors_gene)):
            str_gene = vectors_gene[str_gene_index]
            expression_vector[gene_index_dico[str_gene]] += (radius - 1 * vector_distance[str_gene_index]) / radius
        list_expression_vec.append(expression_vector)
    count_matrix = np.array(list_expression_vec)

    return count_matrix

def get_dict_proba_edge_in_situ(selected_genes, count_matrix,
                                distance="pearson",
                                ):
    """
    compute the colon-expression correlation matrix from the count_matrix
    :param count_matrix: cell by gene matrix
    :type count_matrix: np.array
    :param distance:  choose in ["pearson", "spearman"] default is pearson
    :type distance: str
    :return: a dictionary of dictionary corelation between genes dict[gene_source][gene_target] = correlation
    :rtype: dict
    """
    import math
    assert distance in ["spearman", "pearson"]
    dico_proba_edge = {}
    for gene_source in range(len(selected_genes)):  # I compute the same thing twice ...
        dico_proba_edge[selected_genes[gene_source]] = {}
    for gene_source in tqdm(range(len(selected_genes))):  # I compute the same thing twice ...
        #print(gene_source)
        for gene_target in range(gene_source, len(selected_genes)):
            exp_gene_source = count_matrix[:, gene_source]
            exp_gene_target = count_matrix[:, gene_target]
            if distance == "pearson":
                corr = scipy.stats.pearsonr(exp_gene_source, exp_gene_target)[0]
            elif distance == "spearman":
                corr = scipy.stats.spearmanr(exp_gene_source, exp_gene_target)[0]
            else:
                raise Exception(f'distance {distance} not implemented')
            if math.isnan(corr):
                corr = -1
            dico_proba_edge[selected_genes[gene_source]][selected_genes[gene_target]] = corr
            dico_proba_edge[selected_genes[gene_target]][selected_genes[gene_source]] = corr
    return dico_proba_edge



def compute_coexpression_matrix(list_df_spots_path,
                                list_gene,
                                n_neighbors=40,
                                radius=10,
                                remove_self_node=True,
                                dict_scale={"x": 1, "y": 1, "z": 1},
                                sampling = True,
                                sampling_size = 1000,
                                z_column_name="z",
                                y_column_name="y",
                                x_column_name="x"
                                ):
    """
    compute the colon-expression vector for each gene taking into account all the spots in the list_df_spots_path
    :param list_df_spots_path: list of path to the spots dataframe
    :type list_df_spots_path: list
    :param list_gene: list of gene
    :type list_gene: list
    :param n_neighbors: maximum number of neighbors default is D
    :type n_neighbors: int
    :param radius: maximum radius of neighbors. It should be set proportionnaly to expected cell size, default is radius =  mean_cell_diameter / 2
    :type radius: int
    :param remove_self_node: remove the self edge in the graphe
    :type remove_self_node: bool
    :param dict_scale: dictionnary with the scaling factor for the x,y,z coordinates
    :type dict_scale: dict
    :param sampling: if True sample the count matrix to compute the correlation
    :type sampling: bool
    :param sampling_size: number of vectors to sample from the count matrix to compute correlation
    :type sampling_size: int
    """
    list_of_count_matrix = []
    for df_spots_path in tqdm(list_df_spots_path):
        df_spots = pandas.read_csv(df_spots_path)
        count_matrix = count_matrix_in_situ_from_knn(
                                    list_gene = list_gene,
                                  df_spots_label = df_spots,
                                  n_neighbors=n_neighbors,
                                  radius=radius,
                                  remove_self_node = remove_self_node,
                                  dict_scale = dict_scale,
                                    sampling=sampling,
                                    sampling_size=sampling_size,
                                    z_column_name=z_column_name,
                                    y_column_name=y_column_name,
                                    x_column_name=x_column_name)
        list_of_count_matrix.append(count_matrix)
    ## concatenate the count_matrix
    list_of_count_matrix = np.concatenate(list_of_count_matrix, axis=0)
    if sampling:
        if len(list_of_count_matrix) > sampling_size:
            print("count_matrix.shape", list_of_count_matrix.shape)
            print(f"sampling {sampling} vectors")
            list_of_count_matrix = list_of_count_matrix[np.random.choice(list_of_count_matrix.shape[0], sampling_size, replace=False), :]
            print("count_matrix.shape", list_of_count_matrix.shape)
    dict_co_expression = get_dict_proba_edge_in_situ(list_gene, list_of_count_matrix,
                                distance="pearson",)
    corr_matrix = []
    dict_vect_gene = {}
    for gene0 in dict_co_expression:
        list_corr_gene0 = []
        for gene1 in dict_co_expression:
            list_corr_gene0.append(dict_co_expression[gene0][gene1])
        corr_matrix.append(list_corr_gene0)
        dict_vect_gene[gene0] = list_corr_gene0
    return corr_matrix, dict_vect_gene, dict_co_expression


#%%
if __name__ == "__main__":



    ### transform this script into a command line tool with argparse
    import argparse
    from rna2seg.dataset_zarr.gene2color.gene2color import (
        get_gene_pca_vector,
        get_concatenated_gene_vector,
        get_one_hot_encoded_vector,
    )
    from rna2seg.dataset_zarr.gene2color.utils import get_gene_random_vector

    parser = argparse.ArgumentParser()


    parser.add_argument("--path_to_save", type=str,
                        default="/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/2209_dict_color")
    parser.add_argument('--path_to_cache', type=str,
                        default="/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr/.sopa_cache/rna_seg_1200",
                        help='Path to the nuclei tif file')
    parser.add_argument("--list_df_spots_path", type = str, nargs = "+",  default=None,
                        help="list of path to the spots dataframe")
    parser.add_argument("--list_gene", type = str, nargs = "+",
                        help="list of gene")
    parser.add_argument("--n_neighbors", type = int, default = 40,
                        help="number of neighbors to commpute correlation")
    parser.add_argument("--radius", type = int, default = 10,
                        help="radius of the neighbors to compute correlation")
    parser.add_argument("--remove_self_node", type = bool, default = True)
    parser.add_argument("--dict_scale",
                        #type = dict,
                        default = {"x": 1, "y": 1, "z": 1})
    parser.add_argument("--sampling", type = bool, default = True)
    parser.add_argument("--sampling_size", type = int, default = 200000,
                        help="number of vectors to sample from the count matrix to compute correlation")
    parser.add_argument("--nb_pc", type = int, default = 3,
                        help = "number of principal component to compute gene color, it set the dimention of the gene vector")
    parser.add_argument("--nb_pc_random", type = int, default = 3,
                        help = "The dimention of the random gene vector")
    parser.add_argument("--compute_corr_matrix", type=int, default=1)
    parser.add_argument("--corr_matrix_path", type=str, default=None,
                        help="Used only if compute_corr_matrix is set to 0.")

    parser.add_argument('--nb_crop_sample', type=int, default=60,
                        help='Path to the spots CSV file')

    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')
    args = parser.parse_args()

    Path(args.path_to_save).mkdir(exist_ok=True)


    if args.list_df_spots_path is None:
        list_df_spots_path=[]
        list_path_index = list(Path(args.path_to_cache).glob("*"))
        assert len(list_path_index) > 0, f"no csv found in {args.path_to_cache}"

        # random shuffle the list of index
        np.random.shuffle(list_path_index)
        for path_index in tqdm(list_path_index):
            path_df = path_index  / RNA2segFiles.TRANSCRIPTS_FILE
            if not path_df.exists():
                continue
            df = pd.read_csv(path_df)
            if len(df) > args.n_neighbors:
                list_df_spots_path.append(path_df)
            if len(list_df_spots_path) == args.nb_crop_sample:
                break
            print(f"len(list_df_spots_path) {len(list_df_spots_path)}")


    args.list_df_spots_path = list_df_spots_path
    list_gene = []
    for i in range(len(args.list_df_spots_path)):
        df = pd.read_csv(args.list_df_spots_path[i])
        list_gene += list(df.gene.unique())
    args.list_gene = np.unique(list_gene)


    ## save parameters in the folder
    with open(Path(args.path_to_save) / "script_parameter.txt", "w") as f:
        for k, v in args.__dict__.items():
            f.write(f"{k} : {v}\n")

    if args.compute_corr_matrix:
        corr_matrix, dict_vect_gene, dict_co_expression = compute_coexpression_matrix(
                            list_df_spots_path=args.list_df_spots_path,
                            list_gene=args.list_gene,
                                n_neighbors=args.n_neighbors,
                                radius=args.radius,
                                remove_self_node=args.remove_self_node,
                                dict_scale=args.dict_scale,
                                sampling=args.sampling,
                                sampling_size=args.sampling_size
        )

        ### plot coorr matrix
        list_gene = list(dict_co_expression.keys())
        # plotting the heatmap for correlation
        #ax = sns.heatmap(corr_matrix, xticklabels=list_gene, yticklabels=list_gene, )

        #plt.show()

        np.save(Path(args.path_to_save)/ "corr_matrix", corr_matrix)
    else: 
        corr_matrix = np.load(args.corr_matrix_path)


    # Create gene2color vector
    gene2color = get_gene_pca_vector(coexpression_matrix = corr_matrix,
                                        list_gene  = args.list_gene,
                                        nb_pc = args.nb_pc,
                                        scaling_before_pca=True,
                                        normalize_after_pca=True)
    np.save(Path(args.path_to_save) / f"gene2color_npc{args.nb_pc}", gene2color)

    gene2color_random = get_gene_random_vector(list_gene=list_gene, nb_pc=args.nb_pc_random)
    np.save(Path(args.path_to_save) / f"gene2color_random{args.nb_pc_random}", gene2color_random)

    gene2color_ohe = get_one_hot_encoded_vector(args.list_gene)
    np.save(Path(args.path_to_save) / f"gene2color_ohe", gene2color_ohe)

    concat_gene2color = get_concatenated_gene_vector(gene2color, gene2color_random)
    np.save(Path(args.path_to_save) / f"gene2color_npc{args.nb_pc}_random{args.nb_pc_random}", concat_gene2color)

