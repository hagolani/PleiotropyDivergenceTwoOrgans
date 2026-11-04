# Processed data

Files read by the analysis scripts, built from the raw simulation output.

| File | Content | Used by |
|---|---|---|
| `allexp.dat` | 10 best occluding individuals of the last generation of each simulation: fitness, EMD, cusp numbers | `fig2_analysis.R` |
| `allcuspsallexp.dat` | cusp numbers of upper and lower teeth, all teeth and generations | `fig2_analysis.R` |
| `alldata_monoTooth.dat` | one row per simulation: final fitness, correlation of cusp change (`corr_ul`) | `fig3_analysis.R` |
| `infocategories_threshol0.1_50pbest_delta1_method2.dat` | one row per simulation: inverse fitness, total acceleration | `fig4_analysis.R` |
| `nohoxs{iniend,minmax}.dat`, `hoxs{iniend,minmax}.dat` | parameter ranges per simulation (A_constitutive + hn; hox_cis) | `fig5_analysis.R` |
| `ori_norm_allevos_propsWFIT_exp{1,2,3}.dat` | proportion of adaptive single mutations per simulation and element | `fig6_analysis.R` |
| `boxfit.dat` | final fitness per setting and mutational regime | `sfig_mutation_analysis.R` |

Fitness is stored as inverse fitness in some files; the scripts convert it.
