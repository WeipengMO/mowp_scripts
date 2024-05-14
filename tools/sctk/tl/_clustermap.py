import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import scanpy as sc


class ClusterMap:
    '''
    A class for creating a cluster map based on a given data matrix.

    Example:
    --------
    >>> clustermap = sctk.tl.ClusterMap(corr_matrix)
    >>> clustermap.get_clusters(threshold=.45)
    >>> clustermap.show_clusters()
    >>> clustermap.plot(vmax=.5, vmin=-.5)
    '''

    def __init__(self, data: pd.DataFrame, method: str = 'ward'):
        '''
        Initialize a ClusterMap object.

        Parameters:
        ----------
        data: pd.DataFrame
            A data matrix to be clustered. Usually a correlation matrix.
        
        method: str
            The linkage method to be used. Default is 'ward'.
        '''
        self.data = data  # corr_matrix

        Z = linkage(self.data, method=method)
        self.Z = Z
    

    def _get_dendrogram(self, threshold=.7, show=True, figsize=(12, 3)):
        '''
        Get the dendrogram of the data matrix.
        '''
        Z = self.Z
        color_threshold = threshold*max(Z[:,2])

        if show:
            plt.figure(figsize=figsize)

        dendrogram_res = dendrogram(
            Z, 
            labels=self.data.index, 
            leaf_rotation=90, 
            color_threshold=color_threshold,
            no_plot=not show)
        self.dendrogram_res = dendrogram_res
        

    def get_clusters(self, threshold=.7, show=True, figsize=(12, 3)):
        '''
        Get the clusters based on the dendrogram.

        Parameters:
        ----------
        threshold: float
            The threshold for cutting the dendrogram. Default is 0.7.
        
        show: bool
            Whether to show the dendrogram. Default is True.
        '''
        self._get_dendrogram(threshold=threshold, show=show, figsize=figsize)
        dendrogram_res = self.dendrogram_res
        
        idx = 0
        last_c = None
        cluster_idx = []
        for c in dendrogram_res['leaves_color_list']:
            if last_c is None:
                last_c = c
            if c == last_c:
                cluster_idx.append(f'C{idx}')
            else:
                idx += 1
                cluster_idx.append(f'C{idx}')
                last_c = c
        clusters= pd.Series(cluster_idx, dendrogram_res['ivl'])

        self.clusters = clusters


    def show_clusters(self):
        cluster = {}
        for k, v in self.clusters.items():
            if v not in cluster:
                cluster[v] = []
            cluster[v].append(k)
        
        return cluster

    
    def _get_cluster_color(self, cluster_pal: list = None):
        clusters = self.clusters
        data = self.data

        clusters_color_unique = clusters.unique()

        if cluster_pal is None:
            if len(clusters_color_unique) <= 20:
                cluster_pal = sc.pl.palettes.default_20
            elif len(clusters_color_unique) <= 28:
                cluster_pal = sc.pl.palettes.default_28
            else:
                cluster_pal = sc.pl.palettes.default_102

        if len(clusters_color_unique) > len(cluster_pal):
            raise ValueError(f'Number of clusters exceeds the number of colors in the palette')
        
        cluster_lut = dict(zip(clusters_color_unique, cluster_pal))
        cluster_color = clusters.map(cluster_lut)
        cluster_color = cluster_color[data.index]

        return cluster_color.values


    def plot(
            self, 
            col_colors: list=None, 
            cmap='RdBu_r',
            vmax=None, vmin=None,
            dendrogram_ratio=0.1,
            cluster_pal: list =None,
            figsize=(5, 5)):
        '''
        Plot the cluster map.

        Parameters:
        ----------
        col_colors: list
            A list of colors for the columns. Default is None.
        
        cmap: str
            The colormap to be used. Default is 'RdBu_r'.
        
        vmax: float
            The maximum value for the colormap. Default is None.
        
        vmin: float
            The minimum value for the colormap. Default is None.
        
        dendrogram_ratio: float
            The ratio of the dendrogram. Default is 0.1.

        cluster_pal: list
            A list of colors for the clusters. Default is None.
        
        figsize: tuple
            The size of the figure. Default is (5, 5).
        '''
        
        Z = self.Z
        cluster_color = self._get_cluster_color(cluster_pal=cluster_pal)

        if col_colors is None:
            col_colors = cluster_color
        else:
            if isinstance(col_colors, (pd.Series, list)):
                if len(col_colors) == len(self.data.index):
                    col_colors = [col_colors, cluster_color]
                else:
                    for i, c in enumerate(col_colors):
                        if len(c) != len(self.data.index):
                            raise ValueError(f'col_colors[{i}] has different length with data')
                    col_colors.append(cluster_color)

        g = sns.clustermap(
            self.data, 
            row_linkage=Z, 
            col_linkage=Z,
            cmap=cmap,
            col_colors=col_colors,
            xticklabels=False, 
            yticklabels=False,
            cbar_kws=dict(orientation='horizontal'),
            dendrogram_ratio=dendrogram_ratio,
            figsize=figsize, vmax=vmax, vmin=vmin)

        g.ax_row_dendrogram.set_visible(False) 
        # g.ax_col_dendrogram.set_visible(False)

        heatmap_bbox = g.ax_heatmap.get_position()
        g.ax_cbar.set_position([heatmap_bbox.xmax-heatmap_bbox.width/5, -heatmap_bbox.ymin, heatmap_bbox.width/5, 0.02])

        return g