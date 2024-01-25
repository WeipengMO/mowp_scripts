import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
import anndata as ad
from loguru import logger
import warnings
import itertools
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns
from matplotlib import pyplot as plt
from collections import defaultdict
import scanpy as sc
warnings.simplefilter("ignore")



def overlap(program1, program2):
    # Calculate the overlap between two programs as the ratio of shared genes to total genes
    shared = len(set(program1).intersection(set(program2))) # Count the number of shared genes
    total = len(set(program1).union(set(program2))) # Count the total number of genes
    return shared / total # Return the overlap ratio


def nmf_k(adata: ad.AnnData, k: int, sample: str, filter_genes: bool = True):
    data = adata.X.T
    genes = adata.var_names
    mask = np.array(list(map(lambda x: not x.startswith(('MT', "RPS", "RPL")), genes)))
    # Perform NMF for a given k and return the top 50 genes for each program
    model = NMF(n_components=k, init="nndsvd", random_state=0) # Create an NMF model
    W = model.fit_transform(data)

    programs = {}
    nmf_scores = {}
    for i in range(k):
        # For each program, get the indices of the top 50 genes based on the coefficients
        # top_genes = np.argsort(W[:, i])[-50: ]
        if filter_genes:
            score = W[:, i][mask]
            top_genes = genes[mask]

            index = np.argsort(score)[-50: ]
            score = score[index]
            top_genes = top_genes[index]

        else:
            top_genes = np.argsort(W[:, i])[-50: ]
            score = W[:, i][top_genes]

        name = f'{sample}_{k}_{i}'
        # Append the gene names to the programs list
        programs[name] = top_genes
        # Append the NMF score to the nmf_scores list
        nmf_scores[name] = score

    programs = pd.DataFrame(programs)
    nmf_scores = pd.DataFrame(nmf_scores)

    return programs, nmf_scores


def get_nmf_programs(adata: ad.AnnData, sample_key: str, sample_names: list, k_start: int = 4, k_end: int = 9, min_cells: int = 50, filter_genes: bool = True):
    programs_list, nmf_scores_list = [], []
    for sample in sample_names:
        _adata = adata[adata.obs[sample_key] == sample]
        logger.info(f'Processing {sample}')

        k_values = np.arange(k_start, k_end + 1)
        if len(_adata.obs_names) < min_cells:
            logger.warning(f'{sample} has less than {min_cells} cells, skipping...')
            continue

        for k in k_values:
            _programs, _nmf_scores = nmf_k(_adata, k, sample, filter_genes=filter_genes)
            programs_list.append(_programs)
            nmf_scores_list.append(_nmf_scores)
    
    programs = pd.concat(programs_list, axis=1)
    nmf_scores = pd.concat(nmf_scores_list, axis=1)
    
    return programs, nmf_scores


# def get_recurrent_programs(programs, sample_names, k_start: int = 4, k_end: int = 9):
#     # Check the first criterion: recurs within the tumor
#     # Compare each program with the other programs for the same sample and different k values
#     # If the overlap is at least 0.7, then the program is recurrent

#     recurrent_programs = set()
#     k_values = np.arange(k_start, k_end + 1)
#     for sample in sample_names:
#         logger.info(f'Processing {sample}')
#         for k1 in k_values:
#             for k2 in range(k1 + 1, k_values[-1]+1):
#                 for i in range(k1):
#                     for j in range(k2):
#                         if f'{sample}_{k1}_{i}' in recurrent_programs and f'{sample}_{k2}_{j}' in recurrent_programs:
#                             continue
#                         program1 = programs[f'{sample}_{k1}_{i}']
#                         program2 = programs[f'{sample}_{k2}_{j}']
#                         if overlap(program1, program2) >= 0.7:
#                             recurrent_programs.add(f'{sample}_{k1}_{i}')
#                             recurrent_programs.add(f'{sample}_{k2}_{j}')
    
#     recurrent_programs = sorted(recurrent_programs)
#     return recurrent_programs


def get_recurrent_programs(programs, sample_names, k_start: int = 4, k_end: int = 9):
    # Check the first criterion: recurs within the tumor
    # Compare each program with the other programs for the same sample and different k values
    # If the overlap is at least 0.7, then the program is recurrent

    recurrent_programs = set()
    k_values = range(k_start, k_end + 1)
    for sample in sample_names:
        # logger.info(f'Processing {sample}')
        for k1, k2 in itertools.combinations(k_values, 2):
            for i, j in itertools.product(range(k1), range(k2)):
                program1 = programs[f'{sample}_{k1}_{i}']
                program2 = programs[f'{sample}_{k2}_{j}']
                if overlap(program1, program2) >= 0.7:
                    recurrent_programs.add(f'{sample}_{k1}_{i}')
                    recurrent_programs.add(f'{sample}_{k2}_{j}')
    
    recurrent_programs = sorted(recurrent_programs)
    return recurrent_programs


def get_robust_programs(programs, sample_names, recurrent_programs):
    # Check the second criterion: recurs across tumors
    # Compare each recurrent program with the programs for the other samples and any k value
    # If the overlap is at least 0.2, then the program is robust

    robust_programs = set()
    robust_score = {}
    programs = programs[recurrent_programs]
    recurrent_programs = set(recurrent_programs)

    # Create a dictionary to store the programs for each sample and k value
    sample_programs = {}
    for sample in sample_names:
        sample_programs[sample] = list(filter(lambda x: x.startswith(sample), recurrent_programs))

    # Iterate over all pairs of samples
    for sample1, sample2 in itertools.combinations(sample_names, 2):
        # Iterate over all pairs of programs from the two samples
        for program1_idx, program2_idx in itertools.product(sample_programs[sample1], sample_programs[sample2]):
            assert sample1 != sample2, 'Sample names are the same'
            program1 = programs[program1_idx]
            program2 = programs[program2_idx]
            overlap_score = overlap(program1, program2)
            if overlap_score >= 0.2:
                robust_programs.add(program1_idx)
                robust_score[program1_idx] = max(robust_score.get(program1_idx, 0), overlap_score)
                robust_programs.add(program2_idx)
                robust_score[program2_idx] = max(robust_score.get(program2_idx, 0), overlap_score)

    robust_programs = sorted(robust_programs)
    return robust_programs, robust_score

                 
def get_nonredundant_programs(programs, sample_names, robust_score):
    # Check the third criterion: non-redundant within the tumor
    # Rank the robust programs by their similarity with programs from other tumors and select them in decreasing order
    # If a program is selected, remove any other program within the same tumor that has more than 0.2 overlap with it
    nonredundant_programs = []
    robust_score_idx = sorted(robust_score, key=robust_score.get, reverse=True)
    for sample in sample_names:
        redundant_programs = set()
        filter_idx = list(filter(lambda x: x.startswith(sample), robust_score_idx))
        for idx1 in filter_idx:  # idx = sample_k_i
            if idx1 in redundant_programs:
                continue
            
            nonredundant_programs.append(idx1)
            for idx2 in filter_idx:
                if idx1 == idx2:
                    continue
                if overlap(programs[idx1], programs[idx2]) >= 0.2:
                    redundant_programs.add(idx2)
    
    return nonredundant_programs


def jaccard_similarity(list1, list2):
    set1, set2 = set(list1), set(list2)
    intersection = len(set1.intersection(set2))
    union = len(set1) + len(set2) - intersection
    return float(intersection) / union


def get_jaccard_matrix(programs, nonredundant_programs):
    # Calculate the Jaccard similarity matrix
    size = len(nonredundant_programs)
    jaccard_matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i, size): 
            program1 = programs[nonredundant_programs[i]]
            program2 = programs[nonredundant_programs[j]]
            similarity = jaccard_similarity(program1, program2)
            jaccard_matrix[i, j] = similarity
            if i != j:
                jaccard_matrix[j, i] = similarity  # use symmetry

    return jaccard_matrix


class MetaProgram():
    def __init__(self, adata: ad.AnnData, sample_key: str, sample_names: list, k_start: int = 4, k_end: int = 9, min_cells: int = 50):
        self.adata = adata
        self.sample_key = sample_key
        self.sample_names = sample_names
        self.k_start = k_start
        self.k_end = k_end
        self.min_cells = min_cells

        self.check_min_cells()
        self.check_log()

    
    def check_min_cells(self):
        sample_names = []
        for sample in self.sample_names:
            _adata = self.adata[self.adata.obs[self.sample_key] == sample]
            if len(_adata.obs_names) < self.min_cells:
                logger.info(f'{sample} has less than {self.min_cells} cells, skipping...')
            else:
                sample_names.append(sample)
        self.sample_names = sample_names

    
    def check_log(self):
        if 'counts' in self.adata.layers.keys():
            log_advise = np.allclose(self.adata.X[:10].sum(), self.adata.layers['counts'][:10].sum())
            if log_advise:
                logger.warning('adata.X is not log-normalized, please normalize X')
        else:
            logger.warning('adata.layers["counts"] does not exist, skip checking log')

        if np.any(self.adata.X.A < 0):
            logger.info('adata.X contains negative values, convert the negative values to zero')
            self.adata.X[self.adata.X < 0] = 0

    
    def filter_genes(self):
        adata = self.adata
        logger.info(f'before filtering: {adata.shape[1]} genes')
        self.adata = adata[:, ~adata.var_names.str.startswith(("RPS", "RPL", 'MT'))]
        logger.info(f'after filtering: {self.adata.shape[1]} genes')


    def use_hvg(self, n_top_genes=3000, batch_key=None):
        adata = self.adata
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
        self.adata = adata[:, adata.var['highly_variable']]


    def get_nmf_programs(self, filter_genes=True):
        programs, nmf_scores = get_nmf_programs(self.adata, self.sample_key, self.sample_names, self.k_start, self.k_end, self.min_cells, filter_genes=True)
        recurrent_programs = get_recurrent_programs(programs, self.sample_names, self.k_start, self.k_end)
        robust_programs, robust_score = get_robust_programs(programs, self.sample_names, recurrent_programs)
        nonredundant_programs = get_nonredundant_programs(programs, self.sample_names, robust_score)
        
        self.programs = programs
        self.nmf_scores = nmf_scores
        self.recurrent_programs = recurrent_programs
        self.robust_programs = robust_programs
        self.robust_score = robust_score
        self.nonredundant_programs = nonredundant_programs
    

    def get_meta_programs(self, threshold=0.7, show_dendrogram=True, n_genes=30):
        jaccard_matrix = get_jaccard_matrix(self.programs, self.nonredundant_programs)
        jaccard_matrix_df = pd.DataFrame(jaccard_matrix)

        Z = linkage(jaccard_matrix, method='average')
        color_threshold = threshold*max(Z[:,2])  # default=.7
        dendrogram_res = dendrogram(Z, color_threshold=color_threshold, no_plot=not show_dendrogram)


        idx = 0
        last_c = None
        cluster_idx = []
        for c in dendrogram_res['leaves_color_list']:
            if last_c is None:
                last_c = c
            if c == last_c:
                cluster_idx.append(f'MP{idx}')
            else:
                idx += 1
                cluster_idx.append(f'MP{idx}')
                last_c = c
        row_color= pd.Series(cluster_idx, dendrogram_res['leaves'])
        row_color.name = 'MP'

        mp_consist = defaultdict(lambda: set())
        for mp_id, _nmf_program in zip(cluster_idx, dendrogram_res['leaves']):
            mp_consist[mp_id].add(self.nonredundant_programs[_nmf_program])
        self.mp_consist = mp_consist


        col_pal = sns.husl_palette(len(self.sample_names), s=.8)
        col_lut = dict(zip(self.sample_names, col_pal))

        # Convert the palette to vectors that will be drawn on the side of the matrix
        patient_array = np.array(list(map(lambda x: '_'.join(x.split('_')[:-2]) , self.nonredundant_programs)))
        col_color = pd.Series(patient_array, index=jaccard_matrix_df.columns).map(col_lut)
        col_color.name = 'Patient'

        meta_programs = defaultdict(lambda: [])
        for idx, cluster_name in row_color.to_dict().items():
            meta_programs[cluster_name].append(self.nonredundant_programs[idx])
        
        # To complete the list to topN genes, the genes with the top NMF scores (in either program) were selected
        meta_programs_df = {}
        for cluster_name in meta_programs:
            _meta_programs = self.programs[meta_programs[cluster_name]]
            _meta_programs = _meta_programs.melt(var_name='program', value_name='gene')
            genes = list(_meta_programs['gene'].value_counts().index[:n_genes])
            meta_programs_df[cluster_name] = genes

        meta_programs_df = pd.DataFrame(meta_programs_df)

        self.meta_programs = meta_programs_df
        self.jaccard_matrix_df = jaccard_matrix_df
        self.Z = Z
        self.col_color = col_color

        row_color
        # set row color
        row_color_unique = row_color.unique()
        if len(row_color_unique) <= 20:
            row_pal = sc.pl.palettes.default_20
        elif len(row_color_unique) <= 28:
            row_pal = sc.pl.palettes.default_28
        else:
            row_pal = sc.pl.palettes.default_102
        row_lut = dict(zip(row_color_unique, row_pal))
        row_color = row_color.map(row_lut)
        self.row_color = row_color
    

    def heatmap(self, figsize=(5, 5), vmax=.5, vmin=.1):
        from ..plotting import Colormaps
        colormaps = Colormaps()
    
        g = sns.clustermap(
        self.jaccard_matrix_df, 
        row_linkage=self.Z, col_linkage=self.Z,
        row_colors=self.row_color, col_colors=self.col_color,
        cmap=colormaps['magma_r'], xticklabels=False, yticklabels=False, figsize=figsize, vmax=vmax, vmin=vmin)

        g.ax_row_dendrogram.set_visible(False) 
        g.ax_col_dendrogram.set_visible(False) 

        space = .01
        bbox = g.ax_col_colors.get_position()
        g.ax_col_colors.set_position([bbox.x0, bbox.y0+space, bbox.width, bbox.height])

        bbox = g.ax_row_colors.get_position()
        g.ax_row_colors.set_position([bbox.x0-space, bbox.y0, bbox.width, bbox.height])


        heatmap_bbox = g.ax_heatmap.get_position()
        x0, _y0, _w, _h = g.cbar_pos
        g.ax_cbar.set_position([heatmap_bbox.x0-space*15, heatmap_bbox.y0+heatmap_bbox.height-_h, _w/1.5, _h])
        g.ax_cbar.set_title('Jaccard\nindex', fontsize=10)


        g.ax_heatmap.set_xlabel('NMF programs', fontsize=10)
        g.ax_heatmap.set_ylabel('NMF programs', fontsize=10)

        return g

    
