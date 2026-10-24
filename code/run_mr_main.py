#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM1 main MR pipeline: regenerate all strand-verified main tables in one command.
Prerequisites: the core raw files in data/raw/ (see supplementary/data_manifest.json).
Usage: python run_mr_main.py [--datadir data/raw]
Outputs (results/):
  MR_main_results_strandverified_2026-09-15.tsv            # 24-trait main table (with MDES90)
  IV_table_S8_representative_SNPs_2026-09-15.tsv           # representative IV details
  allele_strand_audit_table_2026-09-15.tsv             # per-IV sign audit
  MR_main_results_standardized_2026-09-15.tsv  # per-SD scale
  excludable_effect_table_strandverified_2026-09-15.tsv        # MDES80/90
"""
import sys, os, csv, math, collections, argparse
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
from scipy.stats import norm

ROOT = os.path.dirname(os.path.abspath(__file__))
SEED_NOTE = 'random seeds: negative controls live in a separate script (see negative_control_calibration_2026-09-15.tsv)'

def main(datadir):
    DR = datadir
    OUT = os.path.join(ROOT, '..', 'results')
    # ---------- 1. EMMAX (UMD3.1.1 raw + liftover ARS index) ----------
    print('[1/6] EMMAX raw table -> allele/statistics dictionary')
    em_meta = {}
    with open(os.path.join(DR, '0gb5mkm04__File_S1_EMMAX_RI_extracted.csv'), encoding='utf-8', errors='replace') as f:
        r = csv.reader(f); h = next(r)
        i_c = h.index('Chromosome'); i_p = h.index('Position'); i_top = h.index('Top Strand'); i_fwd = h.index('Forward Strand')
        i_be = h.index('Predictor Beta'); i_se = h.index('Predictor Beta SE'); i_maf = h.index('Minor Allele D Frequency')
        n = 0
        for row in r:
            try: c = int(row[i_c]); p = int(row[i_p])
            except ValueError: continue
            if 1 <= c <= 29 and p >= 1:
                em_meta[n] = (c, p, row[i_top].strip('[]').split('/'), float(row[i_be]), float(row[i_se]), float(row[i_maf])); n += 1
    ars = {}
    with open(os.path.join(ROOT, '..', 'data', 'derived', 'emmax_arsucd12.bed')) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            try: cnum = int(p[0][3:])
            except ValueError: continue
            i = int(p[3])
            c0, p0, alleles, be, se, maf = em_meta[i]
            ars[(cnum, int(p[1]) + 1)] = (alleles[0], alleles[1], be, se, maf)
    print(f'   EMMAX ARS pairable SNPs: {len(ars)}')

    # ---------- 2. vmcvdnd3q signal grouping ----------
    print('[2/6] vmcvdnd3q fine-mapping CSV grouping')
    sig_rows = collections.defaultdict(list)
    with open(os.path.join(DR, 'vmcvdnd3q__finemapping_summary_stats_dairy_2025.csv'), encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f):
            try:
                sig_rows[(row['Trait'], row['signal'])].append((
                    row['SNPname'], int(row['Chr']), int(row['Pos']), float(row['MAF']),
                    float(row['effect']), float(row['Pval']), float(row['normedProb'] or 0),
                    float(row['R'] or 0), float(row['HWE_Pval'] or 1)))
            except (ValueError, KeyError): continue
    print(f'   signals: {len(sig_rows)}')

    # ---------- 3. S8 signal set ----------
    print('[3/6] S8 high-confidence signals')
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(DR, 'vmcvdnd3q__Supplemental_Tables.xlsx'), read_only=True, data_only=True)
    ws = wb['Table S8']; it = ws.iter_rows(values_only=True)
    hdr = [str(x) if x is not None else '' for x in next(it)]
    si, sni, pi, ti = hdr.index('signal'), hdr.index('SNP name'), hdr.index('normedProb'), hdr.index('Trait')
    s8 = {}
    for row in it:
        if not row or not row[0]: continue
        sig = str(row[si]).strip()
        try: prob = float(row[pi])
        except (TypeError, ValueError): prob = 0.0
        if sig not in s8 or s8[sig]['prob'] < prob:
            s8[sig] = {'trait': str(row[ti]).strip(), 'lead': str(row[sni]).strip(), 'prob': prob}
    print(f'   S8: {len(s8)}')

    EXPOSURE = {'Milk','Fat','Protein','Fat_Percent','Pro_Percent','Stature','Strength','Body_depth','Dairy_form',
                'Rump_angle','Rump_width','Final_score','Rear_legs(side)','Rear_legs(rear)','Foot_angle','Feet_and_legs',
                'Fore_udder_att','Rear_ud_height','Rear_ud_width','Udder_cleft','Udder_depth','Front_teat_pla','Rear_teat_pla','Teat_length'}

    # ---------- 4. representative SNPs (lead or proxy) ----------
    print('[4/6] representative SNPs + allele-strand orientation')
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    rep = []
    for sig, d in s8.items():
        t = d['trait']
        rows = sig_rows.get((t, sig)) or []
        lead_row = next((x for x in rows if x[0] == d['lead']), None)
        chosen = None; is_proxy = False; r_used = 1.0
        if lead_row and (lead_row[1], lead_row[2]) in ars and lead_row[3] >= 0.01 and lead_row[8] > 1e-6:
            chosen = lead_row
        else:
            cands = [x for x in rows if (x[1], x[2]) in ars and x[3] >= 0.01 and x[8] > 1e-6 and x[7] >= 0.5]
            if cands:
                cands.sort(key=lambda x: (-x[7], -x[6]))
                chosen = cands[0]; is_proxy = True; r_used = chosen[7]
        if chosen is None or t not in EXPOSURE: continue
        rep.append(dict(Trait=t, signal=sig, lead=d['lead'], rep_snp=chosen[0], Chr=chosen[1], Pos=chosen[2],
                        is_proxy='Y' if is_proxy else 'N', R_with_lead=r_used, PCIP_lead=d['prob'],
                        MAF_x=chosen[3], effect_x=chosen[4], Pval_x=chosen[5],
                        em_snp=None, em_beta=None, em_se=None, em_af=None, em_a1=None, em_a2=None))
    # allele-strand orientation (effect=Allele2; EMMAX BETA=A1)
    iv2 = []
    audit = []
    stats = collections.Counter()
    for x in rep:
        e = ars[(x['Chr'], x['Pos'])]
        E1, E2, be, se, af = e
        # lead A1/A2 as orientation reference (proxies inherit the lead; limitation stated)
        lead_row = next((y for y in sig_rows[(x['Trait'], x['signal'])] if y[0] == x['lead']), None)
        a1v, a2v = (lead_row[0] and None, None) if False else (None, None)
        # take lead alleles from sig_rows
        lr = next((y for y in sig_rows[(x['Trait'], x['signal'])] if y[0] == x['lead']), None)
        # sig_rows tuple has no allele columns -> light re-read from the original CSV
        # Allele1/Allele2 exist in the CSV - light re-read:
        iv2.append(dict(x, E1=E1, E2=E2, em_beta=be, em_se=se, em_af=af, em_snp=f'{x["Chr"]}:{x["Pos"]}'))
    # allele-letter re-read (light table: only the needed rows)
    need = {(x['Trait'], x['signal']) for x in iv2}
    alleles_map = {}
    with open(os.path.join(DR, 'vmcvdnd3q__finemapping_summary_stats_dairy_2025.csv'), encoding='utf-8', errors='replace') as f:
        for row in csv.DictReader(f):
            k = (row['Trait'], row['signal'])
            if k in need and row['SNPname'] in {y['lead'] for y in iv2 if (y['Trait'], y['signal']) == k}:
                alleles_map[k] = (row['Allele1'], row['Allele2'])
    ivs = []
    for x in iv2:
        k = (x['Trait'], x['signal'])
        if k not in alleles_map: stats['no_alleles'] += 1; continue
        a1v, a2v = alleles_map[k]
        E1, E2 = x['E1'], x['E2']
        if a2v == E1 or a2v == comp.get(E1): sign = 1
        elif a2v == E2 or a2v == comp.get(E2): sign = -1
        else: stats['ambiguous'] += 1; continue
        stats['direct' if sign == 1 else 'flipped'] += 1
        audit.append(dict(Trait=x['Trait'], signal=x['signal'], rep_snp=x['rep_snp'], A1_e=E1, A2_e=E2, A1_v=a1v, A2_v=a2v, sign=sign))
        zx = abs(norm.isf(max(min(x['Pval_x'], 1-1e-300), 1e-300)/2))
        ivs.append(dict(Trait=x['Trait'], signal=x['signal'], rep_snp=x['rep_snp'], is_proxy=x['is_proxy'],
                        R_with_lead=x['R_with_lead'], beta_x=x['effect_x'], Pval_x=x['Pval_x'], F_x=zx*zx,
                        beta_y=sign*x['em_beta'], se_y=x['em_se'], Pval_y=None, MAF_x=x['MAF_x'],
                        MAF_y=min(x['em_af'], 1-x['em_af'])))
    print(f'   IVs (after strand orientation): {len(ivs)}; orientation counts: {dict(stats)}')
    with open(os.path.join(OUT, 'IV_table_S8_representative_SNPs_2026-09-15.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(ivs[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(ivs)
    with open(os.path.join(OUT, 'allele_strand_audit_table_2026-09-15.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(audit[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(audit)

    # ---------- 5. IVW main table ----------
    print('[5/6] IVW main table (with MDES90)')
    def bh(ps):
        m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); out = [0]*m; prev = 1.0
        for rank in range(m, 0, -1):
            i = order[rank-1]; v = ps[i]*m/rank; prev = min(prev, v); out[i] = prev
        return out
    results = []
    for t in sorted({r['Trait'] for r in ivs}):
        rows = [r for r in ivs if r['Trait'] == t]
        W = [r['beta_x']**2/r['se_y']**2 for r in rows]; WD = [r['beta_y']/r['beta_x'] for r in rows]
        k = len(W); sw = sum(W)
        b = sum(w*d for w, d in zip(W, WD))/sw; se = math.sqrt(1/sw)
        p = 2*norm.sf(abs(b/se))
        Q = sum(w*(d-b)**2 for w, d in zip(W, WD)) if k > 1 else 0.0
        loo = []
        if k >= 3:
            for i in range(k):
                W2 = [w for j, w in enumerate(W) if j != i]; WD2 = [d for j, d in enumerate(WD) if j != i]
                sw2 = sum(W2); loo.append(sum(w*d for w, d in zip(W2, WD2))/sw2)
        results.append(dict(Trait=t, n_IV=k, meanF=round(sum(r['F_x'] for r in rows)/len(rows), 1),
                            n_proxy=sum(1 for r in rows if r['is_proxy'] == 'Y'),
                            IVW_beta=b, IVW_se=se, IVW_P=p, Q=Q, Q_df=k-1,
                            LOO_min=min(loo) if loo else float('nan'), LOO_max=max(loo) if loo else float('nan'),
                            MDES80=2.8016*se, MDES90=3.2415*se))
    qs = bh([r['IVW_P'] for r in results])
    for r, q in zip(results, qs): r['BH_q'] = q
    results.sort(key=lambda r: r['IVW_P'])
    with open(os.path.join(OUT, 'MR_main_results_strandverified_2026-09-15.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(results)
    with open(os.path.join(OUT, 'excludable_effect_table_strandverified_2026-09-15.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['Trait','n_IV','IVW_beta','IVW_se','IVW_P','BH_q','MDES80','MDES90'])
        for r in results:
            w.writerow([r['Trait'], r['n_IV'], r['IVW_beta'], r['IVW_se'], r['IVW_P'], r['BH_q'], r['MDES80'], r['MDES90']])

    # ---------- 6. standardized scale ----------
    print('[6/6] standardized-scale table')
    SD = {'Milk': 1100, 'Fat': 35, 'Protein': 30, 'Fat_Percent': 0.35, 'Pro_Percent': 0.22}
    out = []
    for r in results:
        sd = SD.get(r['Trait'])
        if sd:
            b2 = r['IVW_beta']*sd; se2 = r['IVW_se']*sd; m2 = r['MDES80']*sd
            unit = ('per %d kg' % sd) if 'Percent' not in r['Trait'] else ('per %.2f%%' % sd)
        else:
            b2 = r['IVW_beta']; se2 = r['IVW_se']; m2 = r['MDES80']; unit = 'per 1 unit (1-9 scale)'
        out.append(dict(Trait=r['Trait'], exposure_unit=unit, n_IV=r['n_IV'], IVW_beta_sd=b2, IVW_se_sd=se2,
                        IVW_P=2*norm.sf(abs(b2/se2)), BH_q=r['BH_q'], MDES80_sd=m2))
    out.sort(key=lambda x: x['IVW_P'])
    with open(os.path.join(OUT, 'MR_main_results_standardized_2026-09-15.tsv'), 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(out)
    print('\nDONE: main tables regenerated. BH q<0.05:', sum(1 for r in results if r['BH_q'] < 0.05),
          '; nominal P<0.05:', sum(1 for r in results if r['IVW_P'] < 0.05))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--datadir', default=os.path.join(ROOT, '..', 'data', 'raw'))
    a = ap.parse_args()
    main(a.datadir)
