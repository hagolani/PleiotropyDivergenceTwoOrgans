# Builds Supplementary Tables S1-S5 (one sheet each) from the CSV files that
# the analysis scripts write to output/<script name>/.
# Usage: python analysis/make_tables.py   (after Rscript analysis/run_all.R)
import csv, math
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(ROOT, "output") + os.sep
def rd(f): return list(csv.DictReader(open(A+f)))
S = {"1genome_alt0":"With pleiotropy","2genomes_alt0":"Without pleiotropy","2genomeshigh_alt0":"Without pleiotropy, double mutation",
     "1genome":"With pleiotropy","2genomes":"Without pleiotropy","2genomesHigh":"Without pleiotropy, double mutation"}
SETS = ["With pleiotropy","Without pleiotropy","Without pleiotropy, double mutation"]

F  = Font(name="Arial", size=10); FB = Font(name="Arial", size=10, bold=True)
FT = Font(name="Arial", size=12, bold=True); FN = Font(name="Arial", size=9, italic=True, color="555555")
HF = PatternFill("solid", fgColor="D9E1F2"); thin = Side(style="thin", color="999999")

def pval(p):
    p = float(p)
    return "<0.0001" if p < 1e-4 else round(p, 4)
def pname(p):
    if "_" not in p: return p
    a, b = p.split("_", 1)
    b = b.replace("_", " ").replace(".", " ").replace("  ", " ").strip()
    return f"{a} ({b})"
def num(x, d=3): return round(float(x), d)

def sheet(wb, name, title, blocks, notes, widths):
    ws = wb.create_sheet(name); r = 1
    ws.cell(r, 1, title).font = FT; r += 2
    for btitle, header, rows in blocks:
        ws.cell(r, 1, btitle).font = FB; r += 1
        for j, h in enumerate(header, 1):
            c = ws.cell(r, j, h); c.font = FB; c.fill = HF; c.border = Border(bottom=thin)
            c.alignment = Alignment(wrap_text=True, vertical="center")
        r += 1
        for row in rows:
            for j, v in enumerate(row, 1):
                c = ws.cell(r, j, v); c.font = F
                if isinstance(v, float): c.number_format = "0.0000" if abs(v) < 1 else "0.000"
                if j > 1: c.alignment = Alignment(horizontal="right")
            r += 1
        r += 1
    for n in notes:
        ws.cell(r, 1, n).font = FN; r += 1
    for i, w in enumerate(widths, 1): ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = None

wb = Workbook(); wb.remove(wb.active)

# ---------------- Table S1 ----------------
ic2 = rd("fig2_analysis/fig2_stats_initial_conditions.csv")
ic3 = rd("fig3_analysis/fig3_stats_initial_conditions.csv")
ic4 = rd("fig4_analysis/fig4_stats_initial_conditions.csv")
g = {(S[x["setting"]], x["variable"]): x for x in ic2}
k3 = {x["setting"]: x for x in ic3}; k4 = {x["setting"]: x for x in ic4}
rows = []
for s in SETS:
    c, f = g[(s,"cusp_all")], g[(s,"fitness")]
    rows.append([s, num(c["chi2"]), pval(c["p"]), num(f["chi2"]), pval(f["p"]),
                 num(k3[s]["KW_chi2_corr"]), pval(k3[s]["KW_p_corr"]), num(k4[s]["chi2"]), pval(k4[s]["p"])])
rows2 = [[s, num(k3[s]["F_ic"]), pval(k3[s]["p_ic"])] for s in SETS]
sheet(wb, "Table S1", "S1 Table. Effect of initial conditions within each setting.",
      [("Kruskal-Wallis tests, variable ~ initial condition (10 initial conditions, df = 9)",
        ["Setting","Cusp number (upper + lower): chi²","p","Final fitness: chi²","p",
         "Cusp-change correlation: chi²","p","Total acceleration: chi²","p"], rows),
       ("F-test: lm(fitness ~ cusp-change correlation) vs lm(... + initial condition) (SFig corrIC)",
        ["Setting","F","p"], rows2)],
      ["Cusp number and fitness: one value per simulation (mean of its 10 best individuals); fig2_analysis.R (SFig IC).",
       "Cusp-change correlation and F-test: fig3_analysis.R. Total acceleration: fig4_analysis.R.",
       "Kruskal-Wallis p-values are asymptotic; the earlier version used permutation p-values, so p differs slightly while chi² is identical."],
      [38,16,10,14,10,16,10,16,10])

# ---------------- Table S2 ----------------
t2 = [x for x in rd("fig5_analysis/table_S2.csv") if x["parameter"] not in ("generation","row_num")]
order = {"A_constitutive":0,"hox_cis":1,"hox_trans":2}
t2.sort(key=lambda x: (x["range_type"]!="iniend", order[x["element"]]))
EL = {"A_constitutive":"A_constitutive","hox_cis":"hox_cis","hox_trans":"hox_trans"}
blocks = []
for rt, lab in [("iniend","Genetic divergence (Δ begin-end)"),("minmax","Overall genetic exploration (Δ min-max)")]:
    rows = [[pname(x["parameter"]), EL[x["element"]], x["tooth"],
             num(x["delta_with_without"],4), pval(x["p_with_without"]),
             num(x["delta_with_double"],4), pval(x["p_with_double"]), pval(x["p_without_double"])]
            for x in t2 if x["range_type"]==rt]
    blocks.append((lab, ["Parameter","Genetic element","Tooth",
                         "Δ intercept: with − without pleiotropy","p (Tukey)",
                         "Δ intercept: with − without pleiotropy, double mutation","p (Tukey)",
                         "p (Tukey): without vs double mutation"], rows))
sheet(wb, "Table S2", "S2 Table. Fitness-adjusted differences in parameter ranges between settings.", blocks,
      ["Linear model: range ~ final fitness + setting (common slope); Δ intercept = difference between setting coefficients, i.e. the difference in range at equal fitness. Positive values: larger range with pleiotropy.",
       "Ranges min-max normalized per parameter, jointly across settings and genetic elements. INH excluded; Deg and Ntr have no hox_cis element. Source: fig5_analysis.R (table_S2.csv).",
       "Rows 'generation' and 'row_num' in table_S2.csv are data columns, not parameters, and were removed (see fig5_analysis.R fix)."],
      [40,16,8,16,10,18,10,14])

# ---------------- Table S3 ----------------
sh2 = rd("fig2_analysis/fig2_stats_shapiro.csv"); lv2 = {x["variable"]:x for x in rd("fig2_analysis/fig2_stats_levene.csv")}
sh3 = rd("fig3_analysis/fig3_stats_shapiro.csv"); lv3 = {x["variable"]:x for x in rd("fig3_analysis/fig3_stats_levene.csv")}
sh4 = rd("fig4_analysis/fig4_stats_shapiro.csv"); lv4 = rd("fig4_analysis/fig4_stats_levene.csv")[0]
V = [("cusp_all","Cusp number (upper + lower)",sh2,lv2,"Fig 2"),("fitness","Final fitness",sh2,lv2,"Fig 2"),
     ("emd_up","EMD upper tooth",sh2,lv2,"Fig 2"),("emd_low","EMD lower tooth",sh2,lv2,"Fig 2"),
     ("corr_ul","Cusp-change correlation",sh3,lv3,"Fig 3")]
rows = []
for v, lab, sh, lv, fig in V:
    d = {S.get(x["setting"], x["setting"]): x for x in sh if x["variable"]==v}
    rows.append([lab, fig] + sum([[num(d[s]["W"],4), pval(d[s]["p"])] for s in SETS], []) + [num(lv[v]["F"]), pval(lv[v]["p"])])
d4 = {x["setting"]: x for x in sh4}
rows.append(["Total acceleration","Fig 4"] + sum([[num(d4[s]["W"],4), pval(d4[s]["p"])] for s in SETS], []) + [num(lv4["F"]), pval(lv4["p"])])
hdr = ["Variable","Figure"] + sum([[f"{s}: W", "p"] for s in SETS], []) + ["Levene F","p"]
ma = rd("sfig_mutation_analysis/sfig_mutation_stats_assumptions.csv")
RN = [("alt1","High mutation rate (x5)"),("alt2","Low mutation rate (/5)"),("alt6","Genome structure: all equal"),("alt7","Genome structure: cis equal")]
GS = ["1genome","2genomes","2genomeshigh"]
rows_m = []
for a, lab in RN:
    t = {x["test"]: x for x in ma if x["alt"] == a}
    rows_m.append([lab] + sum([[num(t["Shapiro-Wilk "+g_]["statistic"],4), pval(t["Shapiro-Wilk "+g_]["p"])] for g_ in GS], [])
                  + [num(t["Levene"]["statistic"]), pval(t["Levene"]["p"])])
hdr_m = ["Mutational regime"] + sum([[f"{s}: W","p"] for s in SETS], []) + ["Levene F","p"]
sheet(wb, "Table S3", "S3 Table. Tests of statistical assumptions.",
      [("Baseline mutation rates: Shapiro-Wilk normality per setting and Levene homoscedasticity across settings", hdr, rows),
       ("Final fitness under alternative mutational regimes (SFigs mutrate and mutstr)", hdr_m, rows_m)],
      ["Fig 2 variables: one value per individual (10 best per simulation); Fig 3 and Fig 4: one value per simulation.",
       "Alternative regimes: sfig_mutation_analysis.R (sfig_mutation_stats_assumptions.csv). The baseline fitness row there is identical to the Fig 2 row above.",
       "Cusp number under the alternative regimes is not tested in the current scripts; the earlier version of this table had these values from older scripts."],
      [30,8,12,10,12,10,14,10,10,10])

# ---------------- Table S4 ----------------
wx = rd("fig2_analysis/fig2_stats_wilcoxon.csv")
P = {"1genome_alt0 vs 2genomes_alt0":"With vs without pleiotropy",
     "1genome_alt0 vs 2genomeshigh_alt0":"With pleiotropy vs without pleiotropy, double mutation",
     "2genomes_alt0 vs 2genomeshigh_alt0":"Without pleiotropy vs without pleiotropy, double mutation"}
VN = [("fitness","Occlusion (fitness)"),("emd_up","EMD upper"),("emd_low","EMD lower"),("cusp_all","Cusp number (upper + lower)")]
rows = []
for comp in P:
    row = [P[comp]]
    for v,_ in VN:
        x = [y for y in wx if y["variable"]==v and y["comparison"]==comp][0]
        row += [num(x["Z"]), pval(x["p_adj"])]
    rows.append(row)
b2 = ("Fig 2. Wilcoxon rank-sum tests between settings (BH-adjusted p)",
      ["Comparison"] + sum([[f"{n}: Z","p"] for _,n in VN], []), rows)
cr = rd("fig3_analysis/fig3_stats_correlation.csv")
b3 = ("Fig 3. Spearman correlation between final fitness and cusp-change correlation",
      ["Setting","n","r_s","p"], [[x["setting"], int(x["n"]), num(x["r"],4), pval(x["p"])] for x in cr])
an = rd("fig4_analysis/fig4_stats_ancova.csv"); co = {x[""]:x for x in rd("fig4_analysis/fig4_stats_ancova_coefficients.csv")}
b4 = ("Fig 4. Linear models of total acceleration with final fitness as covariate",
      ["Test","df","F","p"], [[ "Setting effect: acc ~ fitness vs acc ~ fitness + setting", int(an[0]["df"]), num(an[0]["F"]), pval(an[0]["p"])],
                              [ "Fitness × setting interaction (equal slopes)", int(an[1]["df"]), num(an[1]["F"]), pval(an[1]["p"])]])
b4c = ("Fig 4. Coefficients of acc ~ fitness + setting (reference: with pleiotropy)",
       ["Term","Estimate","SE","p"], [[n, num(co[k]["Estimate"]), num(co[k]["Std. Error"]), pval(co[k]["Pr(>|t|)"])]
        for k,n in [("fitness","Final fitness (slope)"),("exp2genomes","Without pleiotropy"),("exp2genomesHigh","Without pleiotropy, double mutation")]])
f6 = rd("fig6_analysis/fig6_stats.csv")
b6 = ("Fig 6 / SFig S4. Proportion of adaptive single mutations",
      ["Comparison","Mean difference","p (t-test)","p (Wilcoxon)"],
      [[x["comparison"], num(x["mean_diff"],4), pval(x["p_t"]), pval(x["p_wilcox"])] for x in f6])
mt = rd("sfig_mutation_analysis/sfig_mutation_stats_tests.csv")
RN2 = [("alt0","Baseline")] + RN
PR = [("1genome vs 2genomes","With vs without pleiotropy"),("1genome vs 2genomeshigh","With pleiotropy vs without pleiotropy, double mutation"),
      ("2genomes vs 2genomeshigh","Without pleiotropy vs without pleiotropy, double mutation")]
rows_t = []
for a, lab in RN2:
    t = {x["test"]: x for x in mt if x["alt"] == a}
    rows_t.append([lab, pval(t["Kruskal-Wallis"]["p"])] + [pval(t[k]["p_adj"]) for k,_ in PR])
b7 = ("SFigs mutrate and mutstr. Final fitness between settings under each mutational regime",
      ["Mutational regime","Kruskal-Wallis p"] + [f"{n}: p (Wilcoxon, BH)" for _,n in PR], rows_t)
sheet(wb, "Table S4", "S4 Table. Statistical results for the main and supplementary figures.", [b2,b3,b4,b4c,b6,b7],
      ["Fig 2: one value per individual (10 best per simulation), fig2_analysis.R. Fig 3: fig3_analysis.R. Fig 4: fig4_analysis.R. Fig 6: fig6_analysis.R.",
       "Fig 6: comparisons within a setting are paired (same simulations); comparisons between settings are Welch t-tests and Wilcoxon rank-sum tests.",
       "Mutational regimes: sfig_mutation_analysis.R; in the alternative structures the total mutational influx is rescaled to the baseline.",
       "Fig 5 results: Table S2. Mean differences in Fig 6 are absolute differences in proportion (0.032 = 3.2 percentage points)."],
      [52,14,12,14,12,14,12,16,12])

rc = rd("sfig_raw_analysis/sfig_raw_correlations.csv")
el = lambda q: "hox_trans" if q == "hn" else ("hox_cis" if q.endswith("_hox") else "A_constitutive")
params5 = []
for x in rc:
    if x["parameter"] not in params5: params5.append(x["parameter"])
params5.sort(key=lambda q: order[el(q)])
dd = {(x["setting"], x["parameter"]): x for x in rc}
rows5 = []
for q in params5:
    row = [pname(q[:-4] if q.endswith("_hox") else q), el(q)]
    for st in SETS:
        x = dd[(st, q)]; row += [num(x["r"]), pval(x["p"]), pval(x["p_adj"])]
    rows5.append(row)
sheet(wb, "Table S5", "S5 Table. Correlation between initial and final parameter values.",
      [("Pearson correlation across simulations between the initial value of each initial condition and the final value of each simulation",
        ["Parameter","Genetic element"] + sum([[f"{st}: r","p","p (BH)"] for st in SETS], []), rows5)],
      ["Initial value: population mean at generation 1 of each initial condition, per setting; final value: population mean at the last generation of each simulation.",
       "BH correction across the 41 parameters within each setting. Source: sfig_raw_analysis.R (sfig_raw_correlations.csv); see SFig raw."],
      [40,16,10,10,10,10,10,10,10,10,10])

os.makedirs(os.path.join(ROOT, "output", "tables"), exist_ok=True)
wb.save(os.path.join(ROOT, "output", "tables", "supplementary_tables.xlsx"))
print("Written: output/tables/supplementary_tables.xlsx")

