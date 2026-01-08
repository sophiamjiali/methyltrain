# ==============================================================================
# Script:           defaults.py
# Purpose:          Defines default configurations for the workflow
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-07
# ==============================================================================

# =====| Default Configurations |===============================================

DEFAULT_CONFIG = {

    "project": "TCGA-KIRP",
    "split": [0.60, 0.20, 0.20],        # Train-Validate-Test split
    "seed": 42,

    "array_type": "Illumina Methylation Epic",  # Desired array type
    "genome_build": "hg38",

    "toggles": {

        "average_duplicates": True,   # Average duplicate samples, or take first

        "quality_control": {
            "sample_qc": True,
            "probe_qc": True
        },

        "preprocessing": {
            "imputation": True,       # Perform imputation, else default to 0
            "winsorization": True    
        },

        "aggregation": {
            "batch_correction": True    # Perform batch correction upon cohort
        }
        
    },

    "downloading": {
        "sample_id": "submitter_id",    # Project sample identifier field
        "file_id": "file_name",         # Project file identifier field
        "metadata": [                   # Metadata fields to fetch
            "lalalalalla"
        ]
    },

    "quality_control": {
        "sample_qc": {
            "missing_threshold": 0.05,
            "outlier_threshold": 0
        },

        "probe_qc": {
            "remove_cross_reactive": True,
            "remove_SNP_associated": True,
            "remove_sex_chromosome": True,
            "missing_threshold": 0.05, 
        }
    },

    "preprocessing": {
        "imputation": "lalalalalalallalalalalalalalalallalalalalalalalal",
        "clip_values": [0.001, 99.99] # fix
    }

        
}