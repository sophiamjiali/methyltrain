# ==============================================================================
# Script:           prepare.py.py
# Purpose:          Exposes an end-to-end workflow wrapper function
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

from typing import Dict, List
from pathlib import Path
import anndata as ad

from ..fs.layout import ProjectLayout, CohortLayout
from .steps import (
    download, 
    clean_data, 
    quality_control, 
    preprocess, 
    aggregate, 
    split,
    load_raw_project,
    load_processed_project,
    save_project
)

# =====| Exposed Functions |====================================================

def prepare_dataset(config: Dict, layout: ProjectLayout) -> ad.AnnData:
    """
    Run the full DNA methylation preprocessing workflow on a given project.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    layout : ProjectLayout
        Object representing a project dataset directory layout.

    Returns
    -------
    ad.AnnData
        The processed dataset.
    """
    
    # Ensure directories exist
    layout.initialize()

    download(config, layout)
    clean_data(layout)
    adata = load_raw_project(config, layout)
    adata = quality_control(adata, config, layout)
    adata = preprocess(adata, config)

    return adata
    

def prepare_cohort(config: Dict, 
                   layout: CohortLayout) -> tuple[ad.AnnData, ad.AnnData, 
                                              ad.AnnData]:
    """
    Aggregate the full DNA methylation preprocessing workflow outputs
    from multiple project(s) into a single cohort and split the dataset
    into training, validation, and testing splits.

    Parameters
    ----------
    config : dict
        Configuration dictionary controlling workflow steps.
    processed_paths : list of Path
        Paths to the processed output files (*.h5ad) for each project.
    cohort_dir : Path
        Directory for cohort data.
    """

    # Load each processed project AnnData object 
    project_adatas = [load_processed_project(path) 
                      for path in layout.project_list]
    
    # Aggregate the projects into a cohort AnnData object
    cohort_adata = aggregate(project_adatas)

    # Split the cohort AnnData object into train-val-test splits
    train_adata, val_adata, test_adata = split(cohort_adata, config)

    return train_adata, val_adata, test_adata
