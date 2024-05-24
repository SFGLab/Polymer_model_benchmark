> 4D Nucleome Hackathon March 18-21, 2024
# Project 6: Polymer model benchmarking

## Repository structure

1. We provide information about running the chosen software for chromatin structure modeling in the [run_sims/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/run_sims) folder.
2. We put the scripts for model comparison in the [analysis/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/analysis) folder. Inside we provide the Jupyter Notebooks that constitue the workflow:

- [make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb): to generate distance matrices from output models of chromatin structure (in the XYZ, PDB or CIF format)
- [process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb): to process the distance matrices, calculate Spearman correlation coefficients between those matrices and plot the results of model comparison and validation

Additionally we provide the Conda environment YAML files (one file for each Team member) in the [yamls/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/yamls) folder.

## Workflow

1. Run software for chromatin structure prediction.
2. Interpolate the generated models to the same number of coordinates (e.g., 214).
3. Calculate distance matrices for the models ([make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb)).
4. Compare pairs of distance matrices ([process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb)).

## Data used during the Hackathon

**region of interest** chr1:178421513-179491193

**4DN Data Portal (https://data.4dnucleome.org/)**
- HiC: 4DNES4AABNEZ, 4DNESNMAAN97
- SPRITE:

**ENCODE (https://www.encodeproject.org/)**
- HiC: ENCSR968KAY
- ChIA-PET: ENCSR184YZV_CTCF_ChIAPET (ENCFF379AWZ.hic), ENCSR764VXA
- ChIP-Seq: ENCSR000DZN_CTCF, ENCSR000DZP_SMC

**Rao et al., 2014**
- HiC: GSE63525 (GM12878)

## Software packages used during the Hackathon
- LoopSage (Korsak & Plewczynski 2024)
- MiChroM (Di Pierro et al., 2016)
- DIMES (Shi & Thirumalai, 2023)
- PHi-C2 (Shinkai et al., 2020)
