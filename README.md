## Chromatin model comparison and validation

### 4DN Nucleome Hackathon 2024

On March 18-21, 2024 at The University of Washington in Seattle, USA, we conducted a project to address two challenges in genomic research: comparison and validation of chromatin models. These challenges seem straightforward, however they constitute a difficult challenge. During the hackathon we develop a workflow for model comparison and validation, in which we convert models to distance matrices and we calculate Spearman correlation coefficients between pairs of matrices to estimate the correlations between the models. We followed the workflow with 5 distinct software packages for chromatin modeling. Our results have been published as a preprint on biorxiv[^1].


### Repository structure

- [run_sims](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/run_sims) - instructions how to run certain software for chromatin structure modeling
- [analysis:](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/analysis)

    - [make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb) - generate distance matrices from models (in the XYZ, PDB or CIF format)
    - [process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb) - process distance matrices, calculate Spearman correlation coefficients and plot the results

- [yamls](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/yamls) - conda environment YAML files (one file per software or Team member)

- [scratch](https://github.com/SFGLab/Polymer_model_benchmark/tree/main/scratch) - playground notebooks that we used during the hackathon


### Workflow

1. Run software for chromatin structure prediction and collect output models (in XYZ, PDB or CIF format). If the output constitues an ensemble of models, calculate and proceed with the average over the ensemble.
2. Interpolate the models to the same number of coordinates (e.g., n=214) and calculate distance matrices for each model ([make_dist_mats.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/make_dist_mats.ipynb)).
3. Compare pairs of distance matrices and calculate Spearman correlation coefficients (SCC) ([process_distance_matrices.ipynb](https://github.com/SFGLab/Polymer_model_benchmark/blob/main/analysis/process_distance_matrices.ipynb)).
4. Compare SCC between pairs of matrices to estimate differences between the models.


### Software packages used during the hackathon
- LoopSage[^2]: https://github.com/SFGLab/LoopSage
- MiChroM[^3]: https://open-michrom.readthedocs.io/en/latest/OpenMiChroM.html
- DIMES[^4]: https://github.com/anyuzx/HIPPS-DIMES
- PHi-C2[^5]: https://github.com/soyashinkai/PHi-C2
- MultiMM[^6]: https://github.com/SFGLab/MultiMM


### Data used during the hackathon

**genomic region of interest:** chr1:178421513-179491193

**4DN Data Portal** (https://data.4dnucleome.org/)
- Hi-C: 4DNES4AABNEZ, 4DNESNMAAN97
- SPRITE: 4DNESI1U7ZW9

**ENCODE** (https://www.encodeproject.org/)
- Hi-C: ENCSR968KAY
- ChIA-PET: ENCSR184YZV_CTCF_ChIAPET (ENCFF379AWZ.hic), ENCSR764VXA
- ChIP-Seq: ENCSR000DZN_CTCF, ENCSR000DZP_SMC

**Rao et al., 2014**
- Hi-C: GSE63525 (for Tier 1 cell line GM12878)


### References

[^1]: Biorxiv: https://doi.org/10.1101/2024.10.02.616241

[^2]: Korsak, S., & Plewczynski, D. (2024). LoopSage: An energy-based Monte Carlo approach for the loop extrusion modeling of chromatin. Methods, 223, 106–117. https://doi.org/10.1016/j.ymeth.2024.01.015

[^3]: Di Pierro, M., Zhang, B., Aiden, E. L., Wolynes, P. G., & Onuchic, J. N. (2016). Transferable model for chromosome architecture. Proceedings of the National Academy of Sciences, 113(43), 12168–12173. https://doi.org/10.1073/pnas.1613607113

[^4]: Shi, G., & Thirumalai, D. (2023). A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nature Communications, 14(1), 1150. https://doi.org/10.1038/s41467-023-36412-4

[^5]: Shinkai, S., Itoga, H., Kyoda, K., & Onami, S. (2022). PHi-C2: Interpreting Hi-C data as the dynamic 3D genome state. Bioinformatics, 38(21), 4984–4986. https://doi.org/10.1093/bioinformatics/btac613

[^6]: Korsak, S., Banecki, K., & Plewczynski, D. (2024). Multiscale molecular modelling of chromatin with multimm: From nucleosomes to the whole genome. Cold Spring Harbor Laboratory. http://dx.doi.org/10.1101/2024.07.26.605260
