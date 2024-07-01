import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import scanpy as sc
from typing import Union, Literal


class ClusterMap:
    '''
    A class for creating a cluster map based on a given data matrix.

    Example:
    --------
    >>> clustermap = sctk.tl.ClusterMap(corr_matrix)
    >>> clustermap.get_clusters(cluster = 'both', threshold=.45)
    >>> clustermap.show_clusters()
    >>> clustermap.plot(vmax=.5, vmin=-.5)
    '''

    def __init__(
            self, 
            data: pd.DataFrame,
            ):
        '''
        Initialize a ClusterMap object.

        Parameters:
        ----------
        data: pd.DataFrame
            A data matrix to be clustered. Usually a correlation matrix.
        '''
        self.data = data  # corr_matrix
        self.col_Z = None
        self.row_Z = None
        self.col_clusters = None
        self.row_clusters = None


    def _get_dendrogram(self, Z, labels, threshold=.7, show=True, figsize=(12, 3), dendrogram_kws=None):
        '''
        Get the dendrogram of the data matrix.
        '''
        color_threshold = threshold*max(Z[:,2])

        if show:
            plt.figure(figsize=figsize)

        dendrogram_res = dendrogram(
            Z, 
            labels=labels, 
            leaf_rotation=90, 
            color_threshold=color_threshold,
            no_plot=not show,
            **(dendrogram_kws if dendrogram_kws is not None else {})
        )
        
        return dendrogram_res
        

    def get_clusters(
            self, 
            cluster: Literal['col', 'row', 'both'],
            method: str = 'ward',
            threshold=.7,
            show=True, 
            figsize=(12, 3),
            dendrogram_kws=None):
        '''
        Get the clusters based on the dendrogram.

        Parameters:
        ----------
        cluster: str
            The type of cluster to be performed. Either 'col', 'row', or 'both'.

        method: str
            The linkage method to be used. Default is 'ward'.

        threshold: float
            The threshold for cutting the dendrogram. Default is 0.7.
        
        show: bool
            Whether to show the dendrogram. Default is True.
        
        figsize: tuple
            The size of the figure. Default is (12, 3).
        
        dendrogram_kws: dict
            Additional keyword arguments for the `scipy dendrogram` function.
        '''
        if  cluster == 'both':
            assert self.data.shape[0] == self.data.shape[1], 'The data matrix must be square for both row and column clustering'
            Z = linkage(self.data, method=method)
            labels = self.data.columns
        elif cluster == 'col':
            Z = linkage(self.data.T, method=method)
            labels = self.data.columns
        elif cluster == 'row':
            Z = linkage(self.data, method=method)
            labels = self.data.index
        else:
            raise ValueError(f'cluster must be either "col", "row", or "both"')
        
        dendrogram_res = self._get_dendrogram(
            Z, labels=labels, threshold=threshold, show=show, figsize=figsize, dendrogram_kws=dendrogram_kws)
        
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

        if cluster == 'col':
            self.col_clusters = clusters
            self.col_Z = Z
        elif cluster == 'row':
            self.row_clusters = clusters
            self.row_Z = Z
        elif cluster == 'both':
            self.col_clusters = clusters
            self.col_Z = Z
            self.row_clusters = clusters
            self.row_Z = Z


    def show_clusters(self, cluster = 'col'):
        if cluster == 'col':
            clusters = self.col_clusters
        elif cluster == 'row':
            clusters = self.row_clusters
        else:
            raise ValueError(f'cluster must be either "col" or "row"')
        
        cluster = {}
        for k, v in clusters.items():
            if v not in cluster:
                cluster[v] = []
            cluster[v].append(k)
        
        return cluster

    
    def _get_cluster_color(self, cluster = 'col', cluster_pal: list = None):
        if cluster == 'col':
            clusters = self.col_clusters
        elif cluster == 'row':
            clusters = self.row_clusters
        else:
            raise ValueError(f'cluster must be either "col" or "row"')
        
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
        if cluster == 'col':
            cluster_color = cluster_color[data.columns]
        elif cluster == 'row':
            cluster_color = cluster_color[data.index]

        return cluster_color.values


    def _handle_colors(self, colors, cluster_type, data_dimension, cluster_pal):
        cluster_color = self._get_cluster_color(
            cluster=cluster_type, cluster_pal=cluster_pal) if getattr(self, f"{cluster_type}_clusters") is not None else None
        
        if colors is False:
            return None
        elif colors is None:
            return cluster_color
        elif isinstance(colors, (pd.Series, list, np.ndarray)):
            if len(colors) == data_dimension:
                if cluster_color is None:
                    return colors
                else:
                    return [colors, cluster_color]
            else:
                for i, c in enumerate(colors):
                    if len(c) != data_dimension:
                        raise ValueError(f'{cluster_type}_colors[{i}] has different length with data')
                if cluster_color is None:
                    return colors
                else:
                    return colors + [cluster_color]
        else:
            return colors


    def plot(
        self, 
        col_colors: Union[list, Literal[False]] = None, 
        row_colors: Union[list, Literal[False]] = None,
        cmap='RdBu_r',
        vmax: float = None, 
        vmin: float = None,
        center: float = None,
        dendrogram_ratio: Union[float, tuple] = 0.1,
        col_cluster_pal: list = None,
        row_cluster_pal: list = None,
        show_row_dendrogram: bool = False,
        show_col_dendrogram: bool = False,
        xticklabels=False,
        yticklabels=False,
        x_rotation: float = None,
        y_rotation: float = None,
        figsize: tuple = (5, 5)
    ):
        '''
        Plot the cluster map.

        Parameters:
        ----------
        col_colors : list, optional
            List of colors for each column, used to color the column labels. if False, no color will be used.
        row_colors : list, optional
            List of colors for each row, used to color the row labels. if False, no color will be used.
        cmap : str, optional
            The colormap to use for the heatmap. Default is 'RdBu_r'.
        vmax : float, optional
            The maximum value to use for the colormap. Default is None.
        vmin : float, optional
            The minimum value to use for the colormap. Default is None.
        dendrogram_ratio : float, optional
            The ratio of the dendrogram height to the heatmap height. Default is 0.1.
        col_cluster_pal : list, optional
            List of colors for the column dendrogram clusters. Default is None.
        row_cluster_pal : list, optional
            List of colors for the row dendrogram clusters. Default is None.
        show_row_dendrogram : bool, optional
            Whether to show the row dendrogram. Default is False.
        show_col_dendrogram : bool, optional
            Whether to show the column dendrogram. Default is False.
        figsize : tuple, optional
            The size of the figure. Default is (5, 5).

        Returns:
        -------
        g : seaborn.ClusterGrid
            The cluster grid object representing the plot.

        '''
        
        # Usage for col_colors
        self.col_colors = self._handle_colors(col_colors, 'col', len(self.data.columns), col_cluster_pal)

        # Usage for row_colors
        self.row_colors = self._handle_colors(row_colors, 'row', len(self.data.index), row_cluster_pal)

        col_cluster = False if self.col_Z is None else True
        row_cluster = False if self.row_Z is None else True

        g = sns.clustermap(
            self.data, 
            row_linkage=self.row_Z, 
            col_linkage=self.col_Z,
            cmap=cmap,
            row_cluster = row_cluster,
            col_cluster = col_cluster,
            col_colors=self.col_colors,
            row_colors=self.row_colors,
            xticklabels=xticklabels, 
            yticklabels=yticklabels,
            cbar_kws=dict(orientation='horizontal'),
            dendrogram_ratio=dendrogram_ratio,
            figsize=figsize, 
            vmax=vmax, 
            vmin=vmin,
            center=center)
        
        if not show_row_dendrogram:
            g.ax_row_dendrogram.set_visible(False)
        
        if not show_col_dendrogram:
            g.ax_col_dendrogram.set_visible(False)

        heatmap_bbox = g.ax_heatmap.get_position()
        g.ax_cbar.set_position([heatmap_bbox.xmax-heatmap_bbox.width/5, -heatmap_bbox.ymin, heatmap_bbox.width/5, 0.02])

        if x_rotation is not None:
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=x_rotation, ha='right')
        if y_rotation is not None:
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=y_rotation, ha='right')

        return g