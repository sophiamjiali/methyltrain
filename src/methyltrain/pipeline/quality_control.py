# ==============================================================================
# Script:           quality_control.py
# Purpose:          Internal quality-control functions and implementation logic
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-08
# ==============================================================================

import anndata as ad
import pandas as pd

from typing import Dict

def sample_qc(adata: ad.AnnData, config: Dict) -> ad.AnnData:
    """
    Performs sample-level quality control upon DNA methylation values presented as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix AnnData object with updated metadata suitable for further probe-level quality control and/or preprocessing.

    Steps performed are toggled and configured in the user-provided 
    configurations, including the following options in order:

    1. Remove high missingness (above the provided threshold)
    2. Remove outliers from distribution (above the provided number of SD)
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's raw DNA methylation data at the 
        CpG matrix level.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's sample-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    


    

    return ad.AnnData()

def probe_qc(adata: ad.AnnData, 
             annotation: pd.DataFrame, 
             config: Dict) -> ad.AnnData:
    """
    Performs probe-level quality control upon DNA methylation values presented as a CpG matrix AnnData object. Returns the quality-controlled CpG matrix AnnData object with updated metadata suitable for further preprocessing.

    Steps performed are toggled and ocnfigured in the user-provided 
    configurations, including the following options in order:

    1. Remove cross-reactive probes
    2. Remove SNP-associated probes
    3. Remove multi-mapped probes
    4. Remove sex chromosome probes
    5. Remove high missingness (above the provided threshold)
    
    Parameters
    ----------
    adata : ad.AnnData
        AnnData object representing a project's sample-level quality-controlled 
        DNA methylation data (or raw if QC was not performed) at the CpG matrix 
        level.
    annotation : pd.DataFrame
        Standardized Illumina annotations provided by the package for the 
        project's array type and genome build.
    config : dict
        Configuration dictionary controlling workflow steps.

    Returns
    -------
    ad.AnnData
        AnnData object representing a project's probe-level quality-controlled 
        DNA methylation data at the CpG matrix level with updated metadata.
    """

    return ad.AnnData()