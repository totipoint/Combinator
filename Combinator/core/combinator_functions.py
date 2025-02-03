import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import rankdata
import itertools
import copy
from sklearn.preprocessing import MinMaxScaler
from kneed import KneeLocator
import os
import shutil
import dask.array as da
import scipy
import seaborn as sns
import datetime

def run_combinator(adata: AnnData, combinator_column_name: str, combinator_signature: str, signature: str,project_dir: str, cluster_column: str, exclude_idents: list | None = None, sparse_mode: bool = False, essential_mode: bool = False, temp_dir_keep: bool = False) -> AnnData:
    #documentation
    """
    The main function of Combinator. 
    This function generates an exhaustive list of combinations from a given list of signature genes and infers a pseudotime ranking per cell that is stored in the input AnnData object.

    :param adata: standard scanpy compatible AnnData object that contains all single cell information
    :param combinator_column_name: the column name within the adata.obs will be where Combinator will store its calculated scores
    :param combinator_signature: is the name of the adata.uns key where the signature will be stored
    :param signature: a list of lists that are separated into modules for a cell type signature
    :param project_dir: the directory where the temporary files are created during processing; we don’t specify a default and recommend that users choose a directory on a hard-drive with sufficient space
    :param essential_mode: Boolean to decide if Combinator should be run with certain genes included in every combination of modules
    :param temp_dir_keep: Boolean to decide if temporary directory of processing intermediates is kept or not; by default, set to False meaning processing interemdiates will be erased on completion of script

    :return: AnnData


    """
    
    
#combinator components

def knee_plot(adata, combinator_column_name = None, nonzero = False):
    #documentation
    """
    A helper function that rank-orders the Combinator scores from a given column to aid in deciding a background cutoff.

    :param adata: standard scanpy compatible AnnData object that contains all single cell information after combinator.run has be called on it
    :param combinator_column_name: the column name within the adata.obs that the contain the scores calculated in a previous run of Combinator that will be filtered
    :param nonzero: Boolean to decide whether knee_plot will be constructed on all Combinator score values (default) or just the nonzero values. sometimes, it can be helpful to just look at the nonzero values if there are many cells that are tied for zero to aid in visually identifying the knee 

    :return: a seaborn scatterplot  

    """

def filter_scores(adata, combinator_column_name = None, threshold = None, filtered_combinator_column_name = None):
    """
    Once users have decided upon a cutoff for a given set of Combinator scores, they can use this function to filter those scores and store them into the AnnData object.


    :param adata: standard scanpy compatible AnnData object that contains all single cell information after combinator.run has be called on it
    :param combinator_column_name: the column name within the adata.obs that the contain the scores calculated in a previous run of Combinator that will be filtered
    :param threshold: float indicating the threshold at which Combinator scores in the specified column will be set to zero of they are below it
    :param filtered_combinator_column_name: name of the new Combinator column name that will contain the filtered scores. We recommend keeping the unfiltered original output of Combinator in case users would like to recalculate thresholds
    :return: AnnData 

    """
    
