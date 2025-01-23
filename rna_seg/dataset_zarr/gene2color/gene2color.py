

#%%
import numpy as np

#import matplotlib.pyplot as plt

from sklearn.decomposition import PCA


######### compute the coexpression matrix

def get_gene_pca_vector(coexpression_matrix, list_gene, nb_pc,
                        scaling_before_pca  = True,
                        normalize_after_pca = True):
    """
    compute the pca of the gene coexpression matrix and return the gene2vect dictionary
    Args:
        coexpression_matrix:
        list_gene: name of gene index
        nb_pc:
        scaling_before_pca: apply z = (x - u) / s before pca  https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
        normalize_after_pca: normalize the pca matrix to [0, 1], useful for visualization

    Returns:
        Dictionary mapping gene names to random vectors.
    """

    from sklearn.preprocessing import StandardScaler
    if scaling_before_pca:
        scaler = StandardScaler()
        coexpression_matrix = scaler.fit_transform(coexpression_matrix)
    ## compute pca of "gene coexpression vector"
    pca = PCA(n_components=nb_pc)
    pca.fit(coexpression_matrix)
    pca_matrix = pca.transform(coexpression_matrix)
    for i in range(pca_matrix.shape[1]):
        ## normalize the pca matrix
        if normalize_after_pca:
            pca_matrix[:,i] = ((pca_matrix[:,i] - np.min(pca_matrix[:,i], axis=0)) /
                          (np.max(pca_matrix[:,i], axis=0) - np.min(pca_matrix[:,i], axis=0)))

    ##initilize gene color dictionary
    gene2vect = {}
    for gene in list_gene:
        gene2vect[gene] = np.zeros(nb_pc)
    for gene_index, gene_name in enumerate(list_gene):
        gene2vect[gene_name] += pca_matrix[gene_index, : nb_pc]
    return gene2vect


def get_one_hot_encoded_vector(list_gene):
    """
    One hot encode gene and return gene2vect dictionary.
    Args:
        list_gene: name of gene index

    Returns:
        Dictionary mapping gene names to one hot encoded vectors.
    """

    unique_genes = sorted(set(list_gene))
    num_unique_genes = len(unique_genes)
    
    gene2vect = {}
    for gene in unique_genes:
        one_hot_vector = np.zeros(num_unique_genes, dtype=int)
        one_hot_vector[unique_genes.index(gene)] = 1
        gene2vect[gene] = one_hot_vector
        
    return gene2vect

def get_concatenated_gene_vector(gene2vect_1, gene2vect_2):
    """
    Concatenate two gene2vect dictionaries.
    
    Args:
        gene2vect_1: First gene2vect dictionary.
        gene2vect_2: Second gene2vect dictionary.
        
    Returns:
        Concatenated gene2vect dictionary.
    """
    common_genes = set(gene2vect_1.keys()) & set(gene2vect_2.keys())
    assert len(common_genes)==len(gene2vect_1.keys())
    assert len(common_genes)==len(gene2vect_2.keys())

    gene2vect = {}
    for gene in common_genes:
        if isinstance(gene2vect_1[gene], (int, float)):
            gene2vect[gene] = np.concatenate(([gene2vect_1[gene]], gene2vect_2[gene]))
        elif isinstance(gene2vect_2[gene], (int, float)) :
            gene2vect[gene] = np.concatenate((gene2vect_1[gene], [gene2vect_2[gene]]))
        elif isinstance(gene2vect_1[gene], np.ndarray) and isinstance(gene2vect_2[gene], np.ndarray):
            gene2vect[gene] = np.concatenate(([gene2vect_1[gene]], [gene2vect_2[gene]]))
        else:
            gene2vect[gene] = np.concatenate((gene2vect_1[gene], gene2vect_2[gene]))
    
    return gene2vect


if __name__ == "__main__":
    dict_dist = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene2vect_dist_mean_sigma_standard.npy', allow_pickle=True).item()
    dict_coexpression = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene_color/gene2color_npc4.npy', allow_pickle=True).item()
    dict_coexpression_npc3 = {}
    for gene, vect in dict_coexpression.items():
        dict_coexpression_npc3[gene] = vect[:3]


    np.save('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene2vect_npc3.npy', dict_coexpression_npc3)

    gene2vect_dist_pca_4 = get_concatenated_gene_vector(dict_dist,dict_coexpression)

    np.save('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene2vect_ist_mean_sigma_standard_pca4.npy', gene2vect_dist_pca_4)


    dict_coexpression = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene_color/gene2vect_npc3.npy', allow_pickle=True).item()
    gene2vect_dist_pca_3 = get_concatenated_gene_vector(dict_dist,dict_coexpression)
    np.save('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/lung/dict_color/gene2vect_ist_mean_sigma_standard_pca3.npy', gene2vect_dist_pca_3)


    #######" melanoma ####


    dict_dist = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/dict_color/gene2vect_dist_mean_sigma_standard.npy', allow_pickle=True).item()
    dict_coexpression = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/dict_color/gene2color_npc4.npy', allow_pickle=True).item()
    gene2vect_dist_pca_4 = get_concatenated_gene_vector(dict_dist,dict_coexpression)
    np.save('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/dict_color/gene2vect_ist_mean_sigma_standard_pca4.npy', gene2vect_dist_pca_4)


    dict_coexpression = np.load('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/dict_color/gene2color_npc3.npy', allow_pickle=True).item()
    gene2vect_dist_pca_3 = get_concatenated_gene_vector(dict_dist,dict_coexpression)
    np.save('/home/tom/share/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/dict_color/gene2vect_ist_mean_sigma_standard_pca3.npy', gene2vect_dist_pca_3)
