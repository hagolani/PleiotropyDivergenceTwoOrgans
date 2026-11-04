# Raw simulation output

Not stored in this repository because of its size; download it from Zenodo
(DOI: XXX) and unpack it here. Two scripts read it directly:

- `analysis/sfig_raw_analysis.R` expects
  `<setting>/alt0/Parameters/parameters_done<replicate>_evo<initial condition>.txt`
  (one row per individual and generation; column `names` ends with the
  generation number).
- `analysis/sfig_S2_analysis.R` expects
  `<setting folder>/Dataevos/data<replicate>_evo<initial condition>_best50.txt` and
  `<setting folder>/Paras/paras<1|2>_done<replicate>_evo<initial condition>_best50.dat`
  (top 50% occluding individuals; paras1 = upper tooth, paras2 = lower tooth),
  with setting folders `1genome/fromstart`, `2genomes/fromstart2g`, `2genomes/highmuta`.

Settings: `1genome` = with pleiotropy, `2genomes` = without pleiotropy,
`2genomeshigh` / `highmuta` = without pleiotropy, double mutation.
Replicates (`done`) and initial conditions (`evo`) are numbered from 1.
