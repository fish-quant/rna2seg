import numpy as np


def get_gene_random_vector(list_gene, nb_pc):
    """
    Return a random gene2vect dictionary.
    Args:
        list_gene: name of gene index.
        nb_pc: size of the random vector.
    Returns:
        Dictionary mapping gene names to random vectors.
    """
    gene2vect = {
        gene: np.random.uniform(0, 1, nb_pc)
        for gene in list_gene
    }
    return gene2vect
