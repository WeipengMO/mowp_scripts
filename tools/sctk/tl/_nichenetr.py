from ..utils import rtools
from ..utils.rtools import rcontext
import rpy2.robjects as ro
from rpy2.robjects.packages import importr


R = ro.r


@rcontext
def _load_prior_knowledge(organism, return_env: bool = False, rEnv=None):
    import yaml
    import os

    with open(os.path.dirname(__file__)+'/nichenet_config.ymal') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    config = config[organism]
    lr_network, ligand_target_matrix, weighted_networks = (
        config['lr_network'], 
        config['ligand_target_matrix'], 
        config['weighted_networks'])

    R(f'''
        lr_network <- readRDS('{lr_network}')
        ligand_target_matrix <- readRDS('{ligand_target_matrix}')
        weighted_networks <- readRDS('{weighted_networks}')
    ''')

    if return_env:
        return rEnv['lr_network'], rEnv['ligand_target_matrix'], rEnv['weighted_networks']
  

@rcontext
def run_nichenet_soj(
        seuratObj,
        organism: str,
        sender_celltypes: list,
        receiver_celltypes: str,
        condition_colname: str,
        condition_oi: str,
        condition_reference: str,
        expression_pct = 0.1,
        geneset: str = 'DE',
        lfc_cutoff: float = 0.25,
        rEnv=None,
        ):
    """Run Nichenet analysis on a Seurat object.

    Parameters
    ----------
    seuratObj : Seurat object
        Seurat object.
    cell_type_key : str
        Key in adata.obs that contains cell type information.
    organism : str
        Organism of the data. Must be either 'human' or 'mouse'.
    sender_celltypes : list
        List of sender cell types.
    receiver_celltypes : str
        Receiver cell type.
    condition_colname : str
        Column in adata.obs that contains condition information.
    condition_oi : str
        Condition of interest.
    condition_reference : str
        Reference condition.
    geneset : str
        Geneset to use. Must be either 'DE', 'up' or 'down'.
    lfc_cutoff : float
        Log fold change cutoff.
    rEnv : dict, optional
        Environment to run the R code in.
    """
    importr('Seurat')
    importr('nichenetr')
    importr('tidyverse')

    if isinstance(seuratObj, ro.methods.RS4) and rtools.r2py(R('class')(seuratObj))[0] == 'Seurat':
        rEnv['seuratObj'] = seuratObj
    else:
        raise ValueError('adata must be either a Seurat object.')
    

    rEnv['sender_celltypes'] = rtools.py2r(sender_celltypes)
    rEnv['receiver_celltypes'] = rtools.py2r(receiver_celltypes)
    rEnv['condition_colname'] = rtools.py2r(condition_colname)
    rEnv['condition_oi'] = rtools.py2r(condition_oi)
    rEnv['condition_reference'] = rtools.py2r(condition_reference)
    rEnv['expression_pct'] = rtools.py2r(expression_pct)
    rEnv['geneset'] = rtools.py2r(geneset)
    rEnv['lfc_cutoff'] = rtools.py2r(lfc_cutoff)

    if organism not in ['human', 'mouse']:
        raise ValueError('organism must be either "human" or "mouse"')
    _load_prior_knowledge(organism, rEnv=rEnv)

    R('''
        nichenet_output <- nichenet_seuratobj_aggregate(
            seurat_obj = seuratObj, 
            sender = sender_celltypes, 
            receiver = receiver_celltypes, 
            condition_colname = condition_colname,
            condition_oi = condition_oi,
            condition_reference = condition_reference,
            expression_pct = expression_pct,
            ligand_target_matrix = ligand_target_matrix,
            lr_network = lr_network,
            weighted_networks = weighted_networks,
            geneset = geneset,
            lfc_cutoff = lfc_cutoff,
        )
    ''')

    return rEnv['nichenet_output']


@rcontext
def predict_ligand_activities(
    geneset_oi: list,
    background_expressed_genes: list,
    potential_ligands: list,
    organism: str,
    n_top: int = 30,
    rEnv=None):
    """Predict ligand activities.

    Parameters
    ----------
    geneset_oi : list
        The gene in receiver for which the expression is possibly affected due to communication with other cells.
    background_expressed_genes : list
        The list of genes that are expressed in the receiver.
    potential_ligands : list
        The list of potential ligands in senders.
    n_top : int
        Number of top ligands to return.
    """

    importr('nichenetr')
    importr('tidyverse')

    rEnv['geneset_oi'] = rtools.py2r(geneset_oi)  # target genes of interest
    rEnv['background_expressed_genes'] = rtools.py2r(background_expressed_genes)
    rEnv['potential_ligands'] = rtools.py2r(potential_ligands)
    rEnv['n_top'] = rtools.py2r(n_top)

    if organism not in ['human', 'mouse']:
        raise ValueError('organism must be either "human" or "mouse"')
    _load_prior_knowledge(organism, rEnv=rEnv)

    R('''
        ligand_activities <- predict_ligand_activities(
            geneset = geneset_oi,
            background_expressed_genes = background_expressed_genes,
            ligand_target_matrix = ligand_target_matrix,
            potential_ligands = potential_ligands)

        (ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
            mutate(rank = rank(desc(aupr_corrected))))

        best_upstream_ligands <- ligand_activities %>% top_n(n_top, aupr_corrected) %>%
            arrange(-aupr_corrected) %>% pull(test_ligand)
    ''')

    # Visualization
    R('''
        vis_results <- list()

        ## Infer target genes of top-ranked ligands
        active_ligand_target_links_df <- best_upstream_ligands %>%
            lapply(get_weighted_ligand_target_links,
                    geneset = geneset_oi,
                    ligand_target_matrix = ligand_target_matrix,
                    n = 200) %>% bind_rows()

        active_ligand_target_links <- prepare_ligand_target_visualization(
            ligand_target_df = active_ligand_target_links_df,
            ligand_target_matrix = ligand_target_matrix,
            cutoff = 0.25)

        order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
        order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

        vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

        vis_results$ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Target genes",
                            color = "purple", legend_title = "Regulatory potential") +
        scale_fill_gradient2(low = "whitesmoke",  high = "purple")
        

        ## Infer and receptors of top-ranked ligands
        ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
            best_upstream_ligands, expressed_receptors,
            lr_network, weighted_networks$lr_sig) 

        vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
            ligand_receptor_links_df,
            best_upstream_ligands,
            order_hclust = "both") 

        vis_results$ligand_receptor_network <- (make_heatmap_ggplot(t(vis_ligand_receptor_network), 
            y_name = "Prioritized ligands", x_name = "Receptors ",  
            color = "mediumvioletred", legend_title = "Prior interaction potential"))
    ''')
    
    return rEnv['ligand_activities'], rEnv['vis_results']
