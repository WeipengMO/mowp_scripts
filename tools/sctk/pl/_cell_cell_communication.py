import pandas as pd
import numpy as np
from typing import Literal, Union
from ..pp.adata_process import get_adata_color
from matplotlib import pyplot as plt
from loguru import logger


def pivot_liana_res(
        liana_res: pd.DataFrame, 
        source: str = 'source', 
        target: str = 'target',
        values: str = 'lr_means',
        mode: Literal['counts', 'weight'] = 'counts') -> pd.DataFrame:
    """
    Pivot the liana_res DataFrame to a table with source and target as index and columns, respectively.

    Parameters
    ----------
    liana_res : pd.DataFrame
        The DataFrame containing the results of the liana analysis.
    source : str, optional
        The column name of the source cell type, by default 'source'.
    target : str, optional
        The column name of the target cell type, by default 'target'.
    values : str, optional
        The column name of the values to be aggregated, by default 'lr_means'.
        Only used when mode='weight'.
    mode : Literal['counts', 'weight'], optional
        The mode of the pivot table, by default 'counts'.
        - 'counts': The number of connections between source and target.
        - 'weight': The mean of the values between source and target.
    """
    if mode == 'counts':
        pivot_table = liana_res.pivot_table(index=source, columns=target, aggfunc='size', fill_value=0)
    elif mode == 'weight':
        pivot_table = liana_res.pivot_table(index=source, columns=target, values=values, aggfunc='mean', fill_value=0)

    return pivot_table


def scale_list(arr, min_val=1, max_val=5, keep_zero=True):
    arr = np.array(arr)
    arr_min = np.min(arr)
    arr_max = np.max(arr)

    scaled_arr = (arr - arr_min) / (arr_max - arr_min) * (max_val - min_val) + min_val
    if keep_zero:
        scaled_arr[arr == 0] = 0  # Set the value of 0 to 0
    
    return scaled_arr


def get_mask_df(
        pivot_table: pd.DataFrame, 
        source_cell_type: Union[list, str] = None, 
        target_cell_type: Union[list, str] = None, 
        mode: Literal['and', 'or'] ='or') -> pd.DataFrame:

    if source_cell_type is None and target_cell_type is None:
        return pivot_table

    if isinstance(source_cell_type, str):
        source_cell_type = [source_cell_type]
    if isinstance(target_cell_type, str):
        target_cell_type = [target_cell_type]

    mask_df = pd.DataFrame(np.zeros_like(pivot_table), index=pivot_table.index, columns=pivot_table.columns, dtype=bool)

    if mode == 'or':
        if source_cell_type is not None:
            mask_df.loc[source_cell_type] = True
        if target_cell_type is not None:
            mask_df.loc[:, target_cell_type] = True
    elif mode == 'and':
        if source_cell_type is not None and target_cell_type is not None:
            mask_df.loc[source_cell_type, target_cell_type] = True

    return pivot_table[mask_df].fillna(0)


def align_dataframe(df):
    # 获取行名和列名
    row_names = set(df.index)
    col_names = set(df.columns)
    
    # 找出行名中缺少的列名
    missing_rows = col_names - row_names
    for row in missing_rows:
        df.loc[row] = [0] * len(df.columns)
    
    # 找出列名中缺少的行名
    missing_cols = row_names - col_names
    for col in missing_cols:
        df[col] = [0] * len(df.index)
    
    # 重新排序行和列
    df = df.sort_index().sort_index(axis=1)
    
    return df


def circle_plot(
        adata,
        liana_res: Union[pd.DataFrame, str],
        pivot_mode: Literal['counts', 'weight'] = 'counts',
        source_key: str = 'source',
        target_key: str = 'target',
        values_key: str = 'lr_means',
        cell_type_key: str = 'cell_type',
        source_cell_type: Union[list, str] = None,
        target_cell_type: Union[list, str] = None,
        mask_mode: Literal['and', 'or'] = 'or',
        figure_size: tuple = (5, 5),
        edge_alpha: float = .5,
        edge_arrow_size: int = 10,
        edge_width_scale: tuple = (1, 5),
        node_alpha: float = 1,
        node_size_scale: tuple = (100, 400),
        node_label_offset: tuple = (0.1, -0.2),
        node_label_size: int = 8,
        node_label_alpha: float = .7,
        ):
    """
    Visualize the cell-cell communication network using a circular plot.

    Parameters
    ----------
    adata : AnnData
        The AnnData object.
    liana_res : Union[pd.DataFrame, str]
        The results of the liana analysis. If a string is provided, it will be used as the key to retrieve the results from adata.uns.
    pivot_mode : Literal['counts', 'weight'], optional
        The mode of the pivot table, by default 'counts'.
        - 'counts': The number of connections between source and target.
        - 'weight': The mean of the values between source and target.
    source_key : str, optional
        The column name of the source cell type, by default 'source'.
    target_key : str, optional
        The column name of the target cell type, by default 'target'.
    values_key : str, optional
        The column name of the values to be aggregated, by default 'lr_means'.
        Only used when pivot_mode='weight'.
    cell_type_key : str, optional
        The column name of the cell type, by default 'cell_type'.
    source_cell_type : Union[list, str], optional
        The source cell type to be included in the plot, by default None.
    target_cell_type : Union[list, str], optional
        The target cell type to be included in the plot, by default None.
    mask_mode : Literal['and', 'or'], optional
        The mode of the mask, by default 'or'.
        - 'or': Include the source or target cell type.
        - 'and': Include the source and target cell type.
    figure_size : tuple, optional
        The size of the figure, by default (5, 5).
    edge_alpha : float, optional
        The transparency of the edges, by default .5.
    edge_arrow_size : int, optional
        The size of the arrow, by default 10.
    edge_width_scale : tuple, optional
        The scale of the edge width, by default (1, 5).
    node_alpha : float, optional
        The transparency of the nodes, by default 1.
    node_size_scale : tuple, optional
        The scale of the node size, by default (100, 400).
    node_label_offset : tuple, optional
        The offset of the node label, by default (0.1, -0.2).
    node_label_size : int, optional
        The size of the node label, by default 8.
    node_label_alpha : float, optional
        The transparency of the node label, by default .7.
    """
    
    try:
        import networkx as nx
    except ImportError:
        logger.error('Please install networkx to use this function.')
        return
    
    if isinstance(liana_res, str):
        if liana_res in adata.uns.keys():
            liana_res = adata.uns[liana_res]
        else:
            raise KeyError(f'{liana_res} not found in adata.uns')
    
    pivot_table = pivot_liana_res(
        liana_res, 
        source=source_key, 
        target=target_key, 
        values=values_key, 
        mode=pivot_mode)
    pivot_table = align_dataframe(pivot_table)

    cell_type_colors = get_adata_color(adata, label=cell_type_key)

    # Mask pivot table
    _pivot_table = get_mask_df(
        pivot_table, 
        source_cell_type=source_cell_type, 
        target_cell_type=target_cell_type, 
        mode=mask_mode)

    G = nx.convert_matrix.from_pandas_adjacency(_pivot_table, create_using=nx.DiGraph())
    pos = nx.circular_layout(G)

    edge_color = [cell_type_colors[cell[0]] for cell in G.edges]
    edge_width = np.asarray([G.edges[e]['weight'] for e in G.edges()])
    edge_width = scale_list(edge_width, max_val=edge_width_scale[1], min_val=edge_width_scale[0])

    node_color = [cell_type_colors[cell] for cell in G.nodes]
    node_size = pivot_table.sum(axis=1).values
    node_size = scale_list(node_size, max_val=node_size_scale[1], min_val=node_size_scale[0], keep_zero=False)


    fig, ax = plt.subplots(figsize=figure_size)

    # Visualize network
    nx.draw_networkx_edges(
        G,
        pos,
        alpha=edge_alpha,
        arrowsize=edge_arrow_size,
        arrowstyle='-|>',
        width=edge_width,
        edge_color=edge_color,
        connectionstyle="arc3,rad=-0.3",
        ax=ax
        )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_color,
        node_size=node_size,
        alpha=node_alpha,
        ax=ax
        )
    label_options = {"ec": "k", "fc": "white", "alpha": node_label_alpha}
    _ = nx.draw_networkx_labels(
        G,
        {k: v + np.array(node_label_offset) for k, v in pos.items()},
        font_size=node_label_size,
        bbox=label_options,
        ax=ax
        )

    ax.set_frame_on(False)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    coeff = 1.2
    ax.set_xlim((xlim[0] * coeff, xlim[1] * coeff))
    ax.set_ylim((ylim[0] * coeff, ylim[1]))
    ax.set_aspect('equal')

    return ax


def ccc_dotplot_p9(
    liana_res, 
    colour,  # colour='lr_means'
    size, # size='cellphone_pvals'
    color_range: tuple = None,  # vmax, vmin (0, .32)
    size_range: tuple = None,  # 0.5, 5
    figsize: tuple = None,  # (7, 4)
    cmap_title: str = None,
    size_title: str = None,
    return_fig: bool = True,
    sort_interaction: Literal['target', 'source'] = None,
    interaction_order: list = None,
    cmap='RdBu_r',
    ):
    '''
    Generate a dot plot visualization of interactions between ligands and receptors.

    Parameters:
    -----------
    liana_res : DataFrame
        The dataframe containing the interaction data.
    colour : str
        The column name in `liana_res` to use for coloring the dots.
    size : str
        The column name in `liana_res` to use for sizing the dots.
    color_range : tuple, optional
        The range of values for the color scale. Default is None.
    size_range : tuple, optional
        The range of values for the size scale. Default is None.
    figsize : tuple, optional
        The size of the figure. Default is None.
    cmap_title : str, optional
        The title for the color scale. Default is None.
    size_title : str, optional
        The title for the size scale. Default is None.
    '''
    try:
        from plotnine import ggplot, geom_point, aes, \
            facet_grid, labs, theme_bw, theme, element_text, element_rect, scale_size_continuous, scale_color_cmap
    except ImportError:
        logger.error('Please install plotnine to use this function.')
        return

    if cmap_title is None:
        cmap_title = colour,
    if size_title is None:
        size_title = size

    if sort_interaction == 'target':
        liana_res = liana_res.copy()
        interaction_cat = sorted(set(liana_res['interaction']), key=lambda x: x.split(' -> ')[1])
        liana_res["interaction"] = pd.Categorical(liana_res["interaction"], categories=interaction_cat)
    elif sort_interaction == 'source':
        liana_res = liana_res.copy()
        interaction_cat = sorted(set(liana_res['interaction']), key=lambda x: x.split(' -> ')[0])
        liana_res["interaction"] = pd.Categorical(liana_res["interaction"], categories=interaction_cat)
    elif interaction_order is not None:
        liana_res = liana_res.copy()
        liana_res["interaction"] = pd.Categorical(liana_res["interaction"], categories=interaction_order)

    p = (ggplot(liana_res, aes(x='target', y='interaction', colour=colour, size=size))
            + geom_point()
            + facet_grid('~source')
            + scale_size_continuous(range=size_range)
            + scale_color_cmap(cmap, limits=color_range)
            + labs(color=cmap_title,
                    size=size_title,
                    y="Interactions (Ligand → Receptor)",
                    x="Target",
                    title="Source")
            + theme_bw(base_size=10)
            + theme(
                # legend_text=element_text(size=14),
                # strip_background=element_rect(fill="white"),
                # strip_text=element_text(size=15, colour="black"),
                # axis_text_y=element_text(size=10, colour="black"),
                # axis_title_y=element_text(colour="#808080", face="bold", size=15),
                axis_text_x=element_text(angle=30, ha='right'),
                figure_size=figsize,
                plot_title=element_text(hjust=0.5, size=10)
                )
        )
    
    if return_fig:
        return p
    
    p.draw()


def ccc_dotplot(
        data: pd.DataFrame,
        source: str = 'source',
        target: str = 'target',
        source_target: str = 'source:target',
        interaction: str = 'interaction',
        source_color: str = 'Set1',
        target_color: str = 'Dark2',
        legend_kws: dict = {'fontsize': 8, 'frameon': False},
        sort_interaction: Union[list, str] = ['receptor_complex', 'ligand_complex'],
        sort_interaction_accending: bool = True,
        figsize: tuple = (4, 8),
        swap_axes: bool = False,
        **kwargs,
    ):
    '''
    Generate a dot plot visualization of interactions between ligands and receptors.

    Parameters:
    -----------
    data : DataFrame
        The dataframe containing the interaction data.
    source : str
        The column name of the source cell type.
    target : str
        The column name of the target cell type.
    source_target : str
        The column name of the source:target pair.
    interaction : str
        The column name of the interaction.
    source_color : str
        The color map for the source cell type.
    target_color : str
        The color map for the target cell type.
    legend_kws : dict
        The keyword arguments for the legend.
    sort_interaction : Union[list, str]
        The column name(s) to sort the interactions. Default is ['receptor_complex', 'ligand_complex'].
    sort_interaction_accending : bool
        Whether to sort the interactions in ascending order. Default is True.
    figsize : tuple
        The size of the figure. Default is (4, 8).
    swap_axes : bool
        Whether to swap the axes. Default is False.
    '''
    try:
        from PyComplexHeatmap import HeatmapAnnotation, DotClustermapPlotter, anno_simple, anno_label
    except ImportError:
        logger.error('Please install PyComplexHeatmap to use this function.')
        return

    if isinstance(sort_interaction, str):
        sort_interaction = [sort_interaction]

    
    if not swap_axes:

        # prepare column annotation
        df_col = data[[source, target, source_target]].drop_duplicates().set_index(source_target)

        col_ha = HeatmapAnnotation(
            Source=anno_simple(df_col[source], cmap=source_color, legend=True, add_text=False, legend_kws=legend_kws), 
            Target=anno_simple(df_col[target], cmap=target_color, legend=True, add_text=False, legend_kws=legend_kws),
            verbose=0,
            label_side='left',
            label_kws={'horizontalalignment':'right', 'fontsize': 9},
            axis=1,
            hgap=.75,
        )

        plt.figure(figsize=figsize)
        dm = DotClustermapPlotter(
            data=data, 
            x=source_target, y=interaction, # hue='target', 
            value='lr_means',c='lr_means',s='cellphone_pvals',
            row_cluster=False,
            col_cluster=False,
            cmap='RdBu_r', vmin=0, vmax=.5,
            top_annotation=col_ha,

            # Interactions (ligand -> receptor)
            show_rownames=True,
            row_names_side='left',
            yticklabels_kws={'labelsize': 6},

            y_order=data.sort_values(sort_interaction, ascending=sort_interaction_accending)['interaction'].drop_duplicates().to_list(),

            verbose=0,
            legend_gap=7,
            # spines=False,
            # grid='minor',
            **kwargs,
        )
    
    else:
        # prepare row annotation
        df_row = data[[source, target, source_target]].drop_duplicates().set_index(source_target)

        row_ha = HeatmapAnnotation(
            Source=anno_simple(df_row[source], cmap=source_color, legend=True, add_text=False, legend_kws=legend_kws), 
            Target=anno_simple(df_row[target], cmap=target_color, legend=True, add_text=False, legend_kws=legend_kws),
            verbose=0,
            label_side='top',
            label_kws={'horizontalalignment':'left', 'fontsize': 8, 'rotation': 60},
            axis=0,
            wgap=1,
        )

        plt.figure(figsize=figsize)
        dm = DotClustermapPlotter(
            data=data, 
            x=interaction, y=source_target, # hue='target', 
            value='lr_means',c='lr_means',s='cellphone_pvals',
            row_cluster=False,
            col_cluster=False,
            cmap='RdBu_r', vmin=0, vmax=.5,
            left_annotation=row_ha,

            # Interactions (ligand -> receptor)
            show_colnames=True,
            col_names_side='bottom',
            xticklabels_kws={'labelsize': 6, 'labelrotation': 30},

            x_order=data.sort_values(sort_interaction, ascending=sort_interaction_accending)['interaction'].drop_duplicates().to_list(),

            verbose=0,
            legend_gap=7,
            # spines=False,
            # grid='minor',
            **kwargs,
        )

    return dm