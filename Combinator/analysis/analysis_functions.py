
from sklearn.neighbors import KernelDensity
from scipy import stats
import seaborn as sns; sns.set_theme(color_codes=True)
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, BSpline, InterpolatedUnivariateSpline
import scanpy as sc
from anndata import AnnData
import numpy as np
import pandas as pd
import multiprocessing
import scipy
import itertools
import copy
import os
import shutil
import h5py
import pickle
from scipy.stats import rankdata
import threading
#import concurrent.futures
import numpy_indexed as npi
from sklearn.preprocessing import MinMaxScaler
from kneed import KneeLocator
import os
import shutil
import dask.array as da
import scipy
#from mpi4py import MPI

#import subprocess
import datetime

import csaps
#kwargs tip from here
#https://stackoverflow.com/questions/1496346/passing-a-list-of-kwargs
def simulate_trajectory(adata: AnnData,cluster_column: str, combinator_column_name: str,cell_cutoff: int = 50, smooth: float = 0.15, cells= None, cmap: str | clr.LinearSegmentedColormap = 'magma', **kwargs):
    
    #documentation

    
    
    """ 
    Simulate trajectories along a 100 point differentiation path. 


    :param adata: standard scanpy compatible AnnData object that contains all single cell information after Combinator.run has be called on it
    :param cluster_column: string specifying the column name that contains the cell grouping information that users would like to use to simulate a trajectory through
    :param combinator_column_name: the column name within the adata.obs that the contain the scores calculated in a previous run of Combinator that will be filtered
    :param cell_cutoff: integer specifying the minimum number of cells that each cluster of specified metadata column must have with a Combinator score > 0
    :param smooth: float used by csaps when creating spline through each cluster 
    :param cmap: matplotlib color palette or custom cmap for visualizing the trajectory on the heatmap
    :param kwargs: additional arguments specific for seaborn and matplotlib
    :return: a seaborn clustermap



    """
    
    
    
    


def combine_tajectories(adata: AnnData,scores: list, combined_trajectory_name: str) -> AnnData:
    #documentation

    """
    Given a list of previously calculated Combinator scores, users can use this function to align them and combine them for aid in downstream applications.

    :param adata: standard scanpy compatible AnnData object that contains all single cell information after Combinator.run has be called on it
    :param scores: a standard python list of the names of the Combinator score columns that are stored in adata.obs that you wish to align and combine
    :param combined_trajectory_name: the column name within the adata.obs that will store the combined trajectory scores
    :return: AnnData



    """



def monte_carlo_pvalue(adata: AnnData,cluster_column: str, combinator_column_name: str, combinator_signature: str, project_dir: str, mc_mode: str,analysis_dir: str,mc_permutations: int = 50, cpu_num: int = 5, essential_mode: bool = False, x_percent: float = 0.2, ):
    #documentation

    """
    Run permutation analysis for specified Combinator columns to calculate p-values for each cell's pseudotime ranking using the Monte Carlo method.

    :param adata: standard scanpy compatible AnnData object that contains all single cell information after combinator.run has be called on it
    :param combinator_column_name: is the column name within the adata.obs that the contain the scores calculated in a previous run of Combinator; can be the filtered or unfiltered scores
    :param project_dir: string indicating full path to directory where the user’s project files are. The monte carlo method works by creating temporary file stores for each different permutation that then get removed as results are collected so as to avoid out of memory errors
    :param mc_mode: string indicating which flavor of the algorithm to apply the monte carlo method on. Available choices are:
        ‘mat_wise’: the full dataset is shuffled such that each cell for gene has a random value from the original dataset
        ‘gene_sample’: for each of the signature genes under study, a normal distribution is calculated and positive values are randomly sampled for each cell for that signature. Combinator is then run on each version of the dataset for the specified number of permutations
    :return: numpy array
    

    """
    



