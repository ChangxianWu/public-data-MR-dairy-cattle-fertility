#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VM1 S17: operating-characteristics simulation of uncorrected allele-strand flips in two-sample MR (2026-10-23).

Background: empirically 280/570 (49.3%) of verifiable lead SNPs required orientation flips between the datasets (before/after audit:
10 of 24 trait estimates changed sign). This script generalizes that empirical hazard into a design tool: under an IVW (first-order
weights) framework, randomly flip the recorded orientation of a fraction f of instruments (unaudited analysis) relative to the same data after auditing, quantifying:
(a) type-I error in the null world (expected calibrated - flips do not create false positives, the honest finding);
(b) in the true-effect world (beta = 2.8016 x SE_pooled, the 80%-power design): power loss, estimate attenuation, sign
error rate, false Q heterogeneity alarms, and 95% CI coverage; (c) agreement of unaudited vs audited
estimate signs (matching the empirical 14/24 = 58%) and of significant signals.

Anchor: the empirical flip rate f = 0.493 (237/481 verifiable leads across the 24 exposure traits); k in {5,10,20}
covers the 5-19 IVs per exposure in this design.

Parameterization (z scale): theta_i ~ N(beta, k); with balanced weights the IVW mean has SE exactly 1, so
z = mean(theta) ~ N(beta, 1) and beta is the noncentrality (beta = 2.8016 iff 80% power). Q = sum((theta-mean)^2)/k
is chi-square(k-1). Calibration self-check: at f=0 the type-I error must be ~0.05 and power ~0.80, otherwise the run aborts.

Randomness: the uniform source is a SHA-256 counter-stream deterministic sampler (seed||counter hash bits),
cryptographically strong and exactly reproducible (seed recorded in the ledger and manuscript); normal transform via NormalDist.inv_cdf.
This is a reproducible scientific simulation, not a cryptographic use.

Usage: python strand_flip_simulation.py
Outputs: supplementary/strand_flip_simulation_S17_2026-10-23.tsv + same-named summary JSON
"""
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

OUT_DIR = Path(__file__).resolve().parent.parent / 'supplementary'
N_SIM = 10000
SEED = 20261023
Z_80 = 2.8016          # MDES80 乘子（与主分析一致；此处即非中心性）
Z_ALPHA = 1.959964     # 双侧 0.05
ND = NormalDist()
Z_QCHI = {             # Q 检验 0.05 临界值（df=k-1）——精确卡方分位数
    5: 9.4877, 10: 16.9190, 20: 30.1435,
}
F_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.493, 0.6]   # 0.493 = VM1 实测暴露子集翻转率
K_GRID = [5, 10, 20]


class DetStream:
    """SHA-256 counter-stream deterministic random bit source: cryptographically strong, exactly reproducible with a fixed seed."""

    def __init__(self, seed, stream_id):
        self._base = f"{seed}:{stream_id}".encode("ascii")
        self._ctr = 0
        self._buf = b""

    def _refill(self):
        self._buf = hashlib.sha256(self._base + str(self._ctr).encode("ascii")).digest()
        self._ctr += 1

    def next_uniform(self):
        while len(self._buf) < 8:
            self._refill()
        raw = self._buf[:8]
        self._buf = self._buf[8:]
        return int.from_bytes(raw, "big") / 2 ** 64

    def next_normal(self):
        return ND.inv_cdf(self.next_uniform())


def simulate(beta, f, k, rng, scenario="independent"):
    """One simulation draw: audited vs unaudited analysis of the same instruments.

    scenario="independent": a fraction f of instruments flip independently (Fisher-Yates draw).
    scenario="cluster": with probability f all k instruments of a trait flip together (proxy SNPs inheriting
    the lead orientation plus a dataset-wide coding mismatch). Under cluster flips the unaudited estimate is the mirrored audited estimate:
    |z| is unchanged when flipped (significance preserved, direction reversed); Q and CI width are unchanged.
    """
    per_sd = math.sqrt(k)
    theta = [rng.next_normal() * per_sd + beta for _ in range(k)]
    z = sum(theta) / k                                   # 审计版 IVW z（SE=1）
    if scenario == "cluster":
        flipped_all = rng.next_uniform() < f             # 整簇一锤定音
        theta_t = [-t for t in theta] if flipped_all else theta
    else:
        idx = list(range(k))
        for i in range(k - 1, 0, -1):                    # Fisher-Yates 确定性洗牌
            j = int(rng.next_uniform() * (i + 1))
            idx[i], idx[j] = idx[j], idx[i]
        flipped = set(idx[:round(f * k)])
        theta_t = [-t if i in flipped else t for i, t in enumerate(theta)]
    est_t = sum(theta_t) / k                             # 未审计版估计（z 同尺度）
    q_t = sum((t - est_t) ** 2 for t in theta_t) / k     # 未审计版 Cochran Q ~ χ²(k-1)
    return z, est_t, q_t


def run_cell(beta, f, k, scenario="independent"):
    rng = DetStream(SEED, f"{scenario}_beta{int(beta * 1000)}_f{int(f * 1000)}_k{k}")
    sig_a = sig_u = agree_sign = agree_sig = wrong_sign_u = cover_u = q_reject_u = 0
    ratios = []
    for _ in range(N_SIM):
        z, est_t, q_t = simulate(beta, f, k, rng, scenario)
        ci_lo, ci_hi = est_t - Z_ALPHA, est_t + Z_ALPHA  # SE=1
        if abs(z) >= Z_ALPHA:
            sig_a += 1
        z_t = est_t
        if abs(z_t) >= Z_ALPHA:
            sig_u += 1
        if (z > 0) == (est_t > 0):
            agree_sign += 1                              # 估计符号一致（对应 S15 口径）
        if abs(z) >= Z_ALPHA and abs(z_t) >= Z_ALPHA and (z > 0) == (z_t > 0):
            agree_sig += 1                               # 两侧均显著且同号（保守口径）
        if beta > 0:
            if abs(z_t) >= Z_ALPHA and z_t < 0:
                wrong_sign_u += 1
            if ci_lo <= beta <= ci_hi:
                cover_u += 1
            ratios.append(est_t / beta)
        if q_t > Z_QCHI[k]:
            q_reject_u += 1
    n = N_SIM
    return {
        "scenario": scenario,
        "beta_pooledSE": beta, "flip_fraction": f,
        "eff_flip_frac": round(f * k) / k if scenario == "independent" else f,
        "k_IV": k,
        "typeI_or_power_audited": sig_a / n,
        "typeI_or_power_unaudited": sig_u / n,
        "est_sign_agreement": agree_sign / n,
        "sig_direction_agreement": agree_sig / n,
        "wrong_direction_rate_unaudited": (wrong_sign_u / n) if beta > 0 else "",
        "mean_ratio_est_true": (sum(ratios) / len(ratios)) if beta > 0 else "",
        "coverage95_unaudited": (cover_u / n) if beta > 0 else "",
        "Q_reject_unaudited": q_reject_u / n,
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for beta in (0.0, Z_80):          # 零效应世界 / 80% 功效真效应世界
        for f in F_GRID:
            for k in K_GRID:
                rows.append(run_cell(beta, f, k, "independent"))
    CLUSTER_F_GRID = [0.1, 0.2, 0.3, 0.493, 0.6]     # 整簇翻转（现实惯例错配模式）
    for beta in (0.0, Z_80):
        for f in CLUSTER_F_GRID:
            for k in K_GRID:
                rows.append(run_cell(beta, f, k, "cluster"))

    # 校准自检：f=0 必须复现名义设计，否则参数化有错，立即中止
    for k in K_GRID:
        t1 = next(r["typeI_or_power_audited"] for r in rows
                  if r["beta_pooledSE"] == 0.0 and r["flip_fraction"] == 0.0 and r["k_IV"] == k)
        pw = next(r["typeI_or_power_audited"] for r in rows
                  if r["beta_pooledSE"] == Z_80 and r["flip_fraction"] == 0.0 and r["k_IV"] == k)
        assert 0.04 <= t1 <= 0.06, f"null calibration broken at k={k}: {t1}"
        assert 0.76 <= pw <= 0.84, f"power calibration broken at k={k}: {pw}"

    cols = list(rows[0].keys())
    lines = ['\t'.join(cols)]
    for r in rows:
        lines.append('\t'.join(
            (f"{v:.4f}" if isinstance(v, float) else str(v)) for v in
            [r[c] for c in cols]))
    (OUT_DIR / "strand_flip_simulation_S17_2026-10-23.tsv").write_text(
        '\n'.join(lines) + '\n', encoding="utf-8")

    def pick(beta, f, k, key, scenario="independent"):
        return next(r[key] for r in rows
                    if r["beta_pooledSE"] == beta and r["flip_fraction"] == f and r["k_IV"] == k
                    and r["scenario"] == scenario)

    summary = {
        "design": {"n_sim": N_SIM, "seed": SEED, "rng": "SHA-256 counter stream",
                   "f_grid": F_GRID, "k_grid": K_GRID, "cluster_f_grid": [0.1, 0.2, 0.3, 0.493, 0.6],
                   "anchor_empirical_flip_rate_24exposures": 0.493,
                   "anchor_empirical_sign_agreement_S15": "14/24 = 0.583 (sign changes 10/24)"},
        "null_world_at_empirical_f": {
            "typeI_unaudited_f0.493_k10": pick(0.0, 0.493, 10, "typeI_or_power_unaudited"),
            "est_sign_agreement_f0.493_k10": pick(0.0, 0.493, 10, "est_sign_agreement"),
            "note": "type-I error stays calibrated; the hazard is effect corruption, not false positives",
        },
        "cluster_world_at_empirical_f": {
            "power_unaudited_cluster_f0.493_k10": pick(Z_80, 0.493, 10, "typeI_or_power_unaudited", "cluster"),
            "wrong_direction_significant_cluster_f0.493_k10":
                pick(Z_80, 0.493, 10, "wrong_direction_rate_unaudited", "cluster"),
            "est_sign_agreement_cluster_f0.493_k10":
                pick(Z_80, 0.493, 10, "est_sign_agreement", "cluster"),
            "Q_reject_cluster_f0.493_k10": pick(Z_80, 0.493, 10, "Q_reject_unaudited", "cluster"),
            "note": ("cluster flips preserve significance and reverse direction: the wrong-direction rate "
                     "among all replicates equals f x power"),
        },
        "effect_world_at_empirical_f": {
            "power_audited_k10": pick(Z_80, 0.0, 10, "typeI_or_power_audited"),
            "power_unaudited_f0.493_k10": pick(Z_80, 0.493, 10, "typeI_or_power_unaudited"),
            "mean_ratio_est_true_k10": pick(Z_80, 0.493, 10, "mean_ratio_est_true"),
            "est_sign_agreement_k10": pick(Z_80, 0.493, 10, "est_sign_agreement"),
            "coverage95_unaudited_k10": pick(Z_80, 0.493, 10, "coverage95_unaudited"),
            "Q_reject_unaudited_k10": pick(Z_80, 0.493, 10, "Q_reject_unaudited"),
            "Q_reject_unaudited_f0.2_k10": pick(Z_80, 0.2, 10, "Q_reject_unaudited"),
        },
    }
    (OUT_DIR / "strand_flip_simulation_S17_2026-10-23.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    print("WROTE tsv + json")
    print(json.dumps(summary["null_world_at_empirical_f"], indent=1))
    print(json.dumps(summary["effect_world_at_empirical_f"], indent=1))


if __name__ == "__main__":
    main()
