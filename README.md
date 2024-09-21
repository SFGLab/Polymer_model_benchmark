> 4D Nucleome Hackathon March 18-21, 2024
# Project 6: Polymer model benchmarking

## Repository structure

1. We provide instructions how to run software for chromatin structure modeling in [run_sims/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/run_sims)
2. We put the scripts for model comparison in [analysis/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/analysis):

- [make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb): to generate distance matrices from models (in the XYZ, PDB or CIF format)
- [process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb): to process distance matrices, calculate Spearman correlation coefficients and plot the results

Additionally we provide the Conda environment YAML files (one file for each Team member) in the [yamls/](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/yamls) folder.

## Workflow

1. Run software for chromatin structure prediction.
2. Interpolate the models to the same number of coordinates (e.g., n=214).
3. Calculate distance matrices ([make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb)).
4. Compare pairs of distance matrices ([process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb)).

## Data used during the Hackathon

**region of interest** chr1:178421513-179491193

**4DN Data Portal (https://data.4dnucleome.org/)**
- HiC: 4DNES4AABNEZ, 4DNESNMAAN97
- SPRITE: 4DNESI1U7ZW9

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
