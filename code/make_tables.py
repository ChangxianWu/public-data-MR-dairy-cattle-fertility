#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Journal Tables 1-6: result TSVs -> manuscript tables (Markdown, DOCX-ready)."""
import sys, os, csv, statistics
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
ROOT = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(ROOT, '..', 'results')
OUT = os.path.join(ROOT, '..', 'tables')
os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

def rd(p):
    return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))

def write_md(name, title, header, rows):
    p = os.path.join(OUT, name + '.md')
    with open(p, 'w', encoding='utf-8') as f:
        f.write(f'## {title}\n\n')
        f.write('| ' + ' | '.join(header) + ' |\n')
        f.write('|' + '---|' * len(header) + '\n')
        for r in rows:
            f.write('| ' + ' | '.join(str(x) for x in r) + ' |\n')
    print('  ', name, f'({len(rows)} rows)')

print('Table 1 data sources')
t1 = [
 ['1', 'Exposure GWAS/fine-mapping (Wang 2025, PMID 40912550)', '50,309 US Holstein bulls', '30 traits; 792 high-confidence fine-mapped signals (BFMAP); 11.29M sequence variants', 'Dryad 10.5061/dryad.vmcvdnd3q (accessed 2026-09-15)', 'ARS-UCD1.2'],
 ['2', 'Outcome GWAS (IndexP60 fertility)', '2,448 US Holstein cows', 'Genome-wide EMMAX statistics; 628,526 SNPs; h2 = 0.165-0.255', 'Dryad 10.5061/dryad.0gb5mkm04 (PMID 36848195; accessed 2026-09-15)', 'UMD3.1.1 (lifted to ARS-UCD1.2)'],
 ['3', 'Direction anchor (DGAT1)', 'Same as #1', 'chr14:603,035-612,781; lead 14:609870 with cross-trait direction chain', 'Contained in #1', 'ARS-UCD1.2'],
 ['4', 'Mammary cis-eQTL (context)', '830 samples', '1,793,931 significant pairs; 7,206 genes; nominal per-chromosome pairs', 'CattleGTEx portal (cattlegtex.kiz.ac.cn; accessed 2026-09-16)', 'ARS-UCD1.2'],
]
write_md('Table1_data_sources', 'Table 1. Data sources', ['#', 'Role', 'Population', 'Content', 'Access', 'Genome build'], t1)

print('Table 2 trait mapping (summary)')
t2 = [
 ['Exposure — milk production', 'Milk, Fat, Protein (3)', 'lead+proxy', '36 IVs'],
 ['Exposure — milk composition', 'Fat%, Protein% (2)', 'lead+proxy', '31 IVs'],
 ['Exposure — body conformation', 'Stature, Strength, Body depth, Dairy form, Rump angle/width, Final score (7)', 'lead+proxy', '72 IVs'],
 ['Exposure — feet/legs', 'Rear legs side/rear, Foot angle, Feet and legs (4)', 'lead+proxy', '16 IVs'],
 ['Exposure — udder morphology', 'Fore/rear udder attachment, udder depth/cleft, teat length/placement (8)', 'lead+proxy', '98 IVs'],
 ['Excluded from exposures', 'Gestation length (reverse-causality risk: pregnancy-related); Net Merit (composite profit index containing fertility weights)', '—', '—'],
 ['Total exposure IV pool', '24 traits', '253 representative SNPs (53 lead + 200 proxy, R>=0.558)', '36+31+72+16+98 = 253'],
]
write_md('Table2_trait_mapping', 'Table 2. Trait role mapping and IV allocation', ['Category', 'Traits', 'IV type', 'IVs'], t2)

print('Table 3 main table (strand-verified)')
mr = rd(os.path.join(W, 'MR_main_results_strandverified_2026-09-15.tsv'))
mr.sort(key=lambda r: r['Trait'])
t3 = []
for r in mr:
    p = float(r['IVW_P'])
    t3.append([r['Trait'], r['n_IV'], r['meanF'], f"{float(r['IVW_beta']):+.4f}", f"{float(r['IVW_se']):.4f}",
               f"{p:.2e}" if p < 0.001 else f"{p:.3f}", f"{float(r['BH_q']):.3f}", f"{float(r['MDES80']):.4f}", f"{float(r['MDES90']):.4f}"])
with open(os.path.join(OUT, 'Table3_main_MR.md'), 'a', encoding='utf-8') as f:
    note = ('Note: beta, SE and MDES are in trait-native units per 1 point of the conformation score '
            '(udder/type/feet traits, 1-9 or point scales) and per 1 kg for milk-yield traits (Milk, Fat, Protein; '
            '1 SD of milk yield ~ 1,100 kg); per-SD standardized values are provided in the source table. '
            'MDES80/90 = minimum detectable effect size at 80%/90% power (2.8016/3.2415 x SE). '
            'Mean F = first-stage strength (F = z^2); no instrument below the conventional weak-instrument threshold (min F = 18.4).')
    f.write('\n' + note + '\n')

print('Table 4 nominal-signal sensitivity matrix')
iv = rd(os.path.join(W, 'IV_table_S8_representative_SNPs_2026-09-15.tsv'))
sign = {}
for r in csv.DictReader(open(os.path.join(W, 'allele_strand_audit_table_2026-09-15.tsv'), encoding='utf-8'), delimiter='\t'):
    sign[(r['Trait'], r['signal'])] = int(r['sign'])
from scipy.stats import norm
NOMINAL = ['Fore_udder_att', 'Front_teat_pla', 'Dairy_form', 'Teat_length']
main = {r['Trait']: r for r in mr}
t4 = []
from scipy.stats import norm as _norm
import math as _math
for t in NOMINAL:
    m = main[t]
    # delta-method weights (from the published table)
    d = next((x for x in rd(os.path.join(W, 'sensitivity_delta_weights_strandverified_2026-09-15.tsv')) if x['Trait'] == t), None)
    # PCIP variants (strand-verified)
    pc = [x for x in rd(os.path.join(W, 'sensitivity_PCIP_weighted_drop_strandverified_2026-09-15.tsv')) if x['Trait'] == t]
    pcip_w = next((x for x in pc if x['variant'] == 'pcip'), None)
    drop = next((x for x in pc if x['variant'] == 'drop_low'), None)
    # lead-only
    lo = next((x for x in rd(os.path.join(W, 'sensitivity_lead_only_2026-09-15.tsv')) if x['Trait'] == t), None)
    # random effects (DL): recomputed from the strand-verified IVs
    sub = [r for r in iv if r['Trait'] == t and (t, r['signal']) in sign]
    w_ = [float(r['beta_y'])/float(r['beta_x']) for r in sub]
    v_ = [(float(r['se_y'])/float(r['beta_x']))**2 for r in sub]
    Wt = [1/x for x in v_]
    wf = sum(a*b for a, b in zip(w_, Wt))/sum(Wt)
    q = sum(b*(a-wf)**2 for a, b in zip(w_, Wt))
    dfree = len(w_) - 1
    c_ = sum(Wt) - sum(x*x for x in Wt)/sum(Wt)
    tau2 = max(0.0, (q-dfree)/c_) if c_ > 0 else 0.0
    denom = sum(x/(1+tau2*x) for x in Wt)
    wre = sum(x*a/(1+tau2*x) for a, x in zip(w_, Wt))/denom
    p_re = 2*_norm.sf(abs(wre/_math.sqrt(1/denom)))
    mb = float(m['IVW_beta'])
    def dirflag(beta):
        return 'Yes' if beta*mb > 0 else 'No'
    t4.append([t, m['n_IV'], f"{mb:+.4f} (P={float(m['IVW_P']):.3f})",
               (f"{float(d['delta_beta']):+.4f} (P={float(d['delta_P']):.3f}) {dirflag(float(d['delta_beta']))}" if d else '—'),
               (f"{float(pcip_w['beta']):+.4f} (P={float(pcip_w['P']):.3f}) {dirflag(float(pcip_w['beta']))}" if pcip_w else '—'),
               (f"{float(drop['beta']):+.4f} (P={float(drop['P']):.3f}) {dirflag(float(drop['beta']))}" if drop else '—'),
               (f"{float(lo['lead_beta']):+.4f} (P={float(lo['lead_P']):.3f}, n={lo['n_IV_lead']}) {dirflag(float(lo['lead_beta']))}" if lo and lo.get('lead_beta') else '—'),
               (f"{wre:+.4f} (P={p_re:.3f}) {dirflag(wre)} [tau2={tau2:.1e}]")])
write_md('Table4_nominal_sensitivity', 'Table 4. Nominal signals: consistency across sensitivity analyses (Yes = direction preserved vs main estimate)',
         ['Trait', 'nIV', 'Main IVW', 'Delta-weight', 'PCIP-weighted', 'Drop low-PCIP', 'Lead-only', 'Random-effects (DL)'], t4)

print('Table 5 excludable magnitudes (from run_mr_main)')
t5 = [[r['Trait'], r['n_IV'], f"{float(r['MDES80']):.4f}", f"{float(r['MDES90']):.4f}"] for r in mr]
write_md('Table5_excludable', 'Table 5. Empirically excludable effect sizes (IndexP60 points per 1-point increase of the exposure score; milk-yield traits in kg per 1 kg — see Table 3 note)',
         ['Trait', 'nIV', 'MDES80', 'MDES90'], t5)

print('Table 6 correlation vs causal contrast')
t6 = [
 ['German Holstein, 34,497 (PMID 38589805)', 'MY vs 7 health traits', 'Own data (not public)', 'Pleiotropy filtering only', 'No', 'rg -0.303; causal ~0'],
 ['Chinese Holstein, >10,000 (PMID 41484547)', '6 milk traits -> genes', 'Own data (not public)', 'Not reported', 'Yes', 'Trait-gene mapping'],
 ['PNAS 2026 (PMID 42258707)', 'Rumen melatonin -> methane', 'Own data', 'Not reported', 'Not reported', 'Causal claim'],
 ['JASBT 2026 (PMID 41923184)', 'ADG -> carcass', 'Own data', 'Not reported', 'Yes', 'Causal claim'],
 ['JASBT 2026 (PMID 41851901)', 'Ileal microbiota -> eating speed/fat', 'Own data', 'Not reported', 'Yes', 'Causal claim'],
 ['Sci China Life Sci 2026 (PMID 41811561)', 'Host-microbe co-regulation -> VFA', 'Own data', 'Not reported', 'Not reported', 'Causal claim'],
 ['Anim Biosci 2026 (PMID 41927046)', 'Placental efficiency -> piglet traits (porcine)', 'Own data', 'Not reported', 'Not reported', 'Causal claim'],
 ['German Holstein GSMR, 235,164 cows (PMID 40000939)', '25 traits: conformation/metabolism/production -> reproduction (156 combos)', 'Own data (not public)', 'Not reported', 'No', '135/156 significant, partially non-linear; within-population (split half-cohorts)'],
 ['Chinese Holstein conformation eQTL-MR, 7,674 (PMID 40383388)', 'Conformation -> causal genes via cis-eQTL MR', 'Own data (not public)', 'Not reported', 'Via eQTL (gene level)', '153 gene-trait MR links'],
 ['This study (VM1)', '24 exposures -> IndexP60', 'Fully public (Dryad x2 + portal)', 'Per-SNP strand audit (280/570 flipped)', 'Proximity enrichment (N=830 eQTL)', 'All BH q>0.05; MDES 0.029'],
]
write_md('Table6_livestock_MR_contrast', 'Table 6. Livestock MR studies: data, strand verification, and findings',
         ['Study', 'Design', 'Data', 'Strand verification', 'Colocalization', 'Key finding'], t6)

print('\nAll 6 tables written ->', OUT)
