import pandas as pd
import numpy as np
import math
from sklearn.preprocessing import MinMaxScaler
from collections import OrderedDict
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import seaborn as sns


class DotPlot:
    def __init__(
            self, 
            value: pd.DataFrame, 
            size: pd.DataFrame, 
        ):
        self.value_df = value.copy()
        self.size_df = size.copy()

    
    @staticmethod
    def _hierarchical_clustering(data: pd.DataFrame):
        Z = linkage(data)
        dendrogram_res = dendrogram(Z, labels=data.index, no_plot=True)
        return dendrogram_res['ivl']
    

    def plot(
            self, 
            reverse_size: bool = True, 
            hspace: float =.2, 
            wspace: float =.2,
            vmax: float = None,
            vmin: float = None,
            max_size: float = None,
            min_size: float = None,
            row_cluster: bool = False, 
            col_cluster: bool = False,
            cmap: str = 'Reds',
            figsize: tuple = (5, 4),
            sizes: tuple = (20, 200),
            cbar_title: str = 'Value',
            size_title: str = 'Size',
            rotation: float = None,
            ):
            '''
            Plot a dotplot using the provided data.

            Parameters:
            -----------
            reverse_size : bool, optional
                Whether to reverse the size values. Default is True.
            hspace : float, optional
                The amount of height reserved for space between subplots, expressed as a fraction of the average axis height. Default is 0.2.
            wspace : float, optional
                The amount of width reserved for space between subplots, expressed as a fraction of the average axis width. Default is 0.2.
            vmax : float, optional
                The maximum value to be shown on the colorbar. Default is None.
            vmin : float, optional
                The minimum value to be shown on the colorbar. Default is None.
            max_size : float, optional
                The maximum size value to be used. Default is None.
            min_size : float, optional
                The minimum size value to be used. Default is None.
            row_cluster : bool, optional
                Whether to cluster the rows. Default is False.
            col_cluster : bool, optional
                Whether to cluster the columns. Default is False.
            cmap : str, optional
                The color map to be used for the dotplot. Default is 'Reds'.
            figsize : tuple, optional
                The figure size. Default is (5, 4).
            sizes : tuple, optional
                The range of sizes for the dots. Default is (20, 200).
            cbar_title : str, optional
                The title for the colorbar. Default is 'Value'.
            size_title : str, optional
                The title for the size legend. Default is 'Size'.
            rotation : float, optional
                The rotation angle for the x-axis labels. Default is None.
            '''

            if row_cluster:
                row_order = self._hierarchical_clustering(self.value_df)
                self.value_df = self.value_df.loc[row_order]
            if col_cluster:
                col_order = self._hierarchical_clustering(self.value_df.T)
                self.value_df = self.value_df[col_order]

            # Prepare the data for plotting
            df = (self.value_df.reset_index().
                  rename(columns={'index': 'col_names'}).
                  melt(id_vars='col_names', var_name='row_names', value_name='Value'))
            size_df = (self.size_df.reset_index().
                       rename(columns={'index': 'col_names'}).
                       melt(id_vars='col_names', var_name='row_names', value_name='Size'))

            if vmax is not None:
                df['Value'] = df['Value'].clip(upper=vmax)
            if vmin is not None:
                df['Value'] = df['Value'].clip(lower=vmin)

            if max_size is not None:
                df['Size'] = df['Size'].clip(upper=max_size)
            if min_size is not None:
                df['Size'] = df['Size'].clip(lower=min_size)

            if reverse_size:
                df['Size'] = -size_df['Size'].apply(lambda x: np.log10(x))


            fig = plt.figure(figsize=figsize)

            # Define the grid layout
            # 通过width_ratios间接控制colorbar的宽度
            gs = fig.add_gridspec(2, 2, width_ratios=[20, 1], height_ratios=[1, 1], hspace=hspace, wspace=wspace)

            # Scatter plot
            ax1 = fig.add_subplot(gs[:, 0])
            g = sns.scatterplot(
                data=df, 
                x='col_names', 
                y='row_names', 
                size='Size', 
                sizes=sizes, 
                hue='Value',
                palette=cmap, 
                edgecolor="w", 
                legend=False,
                ax=ax1
            )

            g.set_xlabel('')
            g.set_ylabel('')
            
            # Add space for x-axis
            xlim = g.get_xlim()
            space = (xlim[1] - xlim[0]) * 0.1
            g.set_xlim(xlim[0] - space, xlim[1] + space)

            if rotation is not None:
                plt.setp(g.xaxis.get_majorticklabels(), rotation=rotation, ha='right')

            # Colorbar
            ax2 = fig.add_subplot(gs[0, 1])
            norm = plt.Normalize(df['Value'].min(), df['Value'].max())
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, cax=ax2, orientation='vertical')
            cbar.ax.tick_params(labelsize=8) 
            ax2.set_title(cbar_title, fontsize=8)

            # Legend
            ax3 = fig.add_subplot(gs[1, 1])
            ax3.axis('off')
            dotsize = np.linspace(math.floor(df['Size'].min()), math.ceil(df['Size'].max()), 4, dtype=int)
            dotsize = np.array(list(OrderedDict.fromkeys(dotsize)))
            scaler = MinMaxScaler(feature_range=sizes)
            scale_dotsize = scaler.fit_transform(dotsize.reshape(-1, 1)).reshape(-1)            
            for size, label in zip(scale_dotsize, dotsize):
                ax3.scatter([], [], c='k', alpha=0.2, s=size, label=f'{label}')
            
            legend = ax3.legend(scatterpoints=1, frameon=False, labelspacing=1, loc='center')
            legend.set_title(size_title, prop={'size': 8})
            ax3.add_artist(legend)

            plt.show()
