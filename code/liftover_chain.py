#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lightweight liftOver (UCSC chain parser, zero dependencies): EMMAX IndexP60 (UMD3.1.1/bosTau6) -> ARS-UCD1.2 (bosTau9).
Usage:
  python liftover_chain.py build_index <chain.gz> <index.pkl>          # pre-parse the index
  python liftover_chain.py convert <index.pkl> <input.bed> <output.bed> <unmapped.bed>   # 0-based BED conversion
Chain format (UCSC): chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
Data lines: size [dt dq], cumulative within the t chain; minus-strand q runs opposite.
Index: {tName: [(tStart, tEnd, qName, qStart, tStrand, qStrand, blocks), ...]} blocks=[(size, dt, dq), ...]
Conversion: for a query point p on the t chain find the covering chain, locate it in blocks, q = qStart + offset; minus-strand target uses the qEnd offset.
"""
import gzip, sys, pickle, bisect, os

def build_index(chain_gz):
    idx = {}
    cur = None
    def flush():
        if cur and cur['blocks']:
            idx.setdefault(cur['tName'], []).append(
                (cur['tStart'], cur['tEnd'], cur['qName'], cur['qStart'], cur['qEnd'], cur['tStrand'], cur['qStrand'], tuple(cur['blocks'])))
    with gzip.open(chain_gz, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('chain'):
                flush()
                p = line.split()
                cur = {'tName': p[2], 'tStart': int(p[5]), 'tEnd': int(p[6]), 'tStrand': p[4],
                       'qName': p[7], 'qStart': int(p[10]), 'qEnd': int(p[11]), 'qStrand': p[11 - 1] if False else p[12 - 1], 'blocks': []}
                # 注意：qStrand 是第 11 个字段（0-based index 11 为 id），实际: p[4]=tStrand? 顺序见 docstring
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                cur['tStrand'] = p[4]; cur['qStrand'] = p[9]
                cur['tStart'] = int(p[5]); cur['tEnd'] = int(p[6])
                cur['qName'] = p[7]; cur['qStart'] = int(p[10]); cur['qEnd'] = int(p[11])
                cur['qStrand'] = p[9]
                cur['blocks'] = []
            else:
                p = line.split()
                size = int(p[0]); dt = int(p[1]) if len(p) > 1 else 0; dq = int(p[2]) if len(p) > 2 else 0
                cur['blocks'].append((size, dt, dq))
    flush()
    # 排序便于二分
    for k in [k for k in idx if not k.endswith('__starts')]:
        idx[k].sort(key=lambda x: x[0])
        idx[k + '__starts'] = [c[0] for c in idx[k]]
        me = 0; maxend = []
        for c in idx[k]:
            me = max(me, c[1]); maxend.append(me)
        idx[k + '__maxend'] = maxend
    return idx

def convert_point(idx, t_name, pos0):
    arr = idx.get(t_name)
    if not arr:
        return None
    starts = idx[t_name + '__starts']
    maxend = idx[t_name + '__maxend']
    # 覆盖 pos0 的 chain = 起点 <= pos0 且 tEnd > pos0；用前缀最大 tEnd 二分最右候选
    hi = bisect.bisect_right(maxend, pos0)  # maxend[hi-1] <= pos0 → 这些 chain 全部 tEnd<=pos0
    # 从 hi 向左找第一个 tStart<=pos0 的 chain（含 gap 跳过）
    j = bisect.bisect_right(starts, pos0) - 1
    while j >= 0:
        tStart, tEnd, qName, qStart, qEnd, tStrand, qStrand, blocks = arr[j]
        if tEnd > pos0 and tStart <= pos0:
            t_acc = tStart
            q_acc = qStart if qStrand == '+' else qEnd
            for (size, dt, dq) in blocks:
                if t_acc <= pos0 < t_acc + size:
                    off = pos0 - t_acc
                    return (qName, q_acc + off) if qStrand == '+' else (qName, q_acc - off)
                t_acc += size + dt
                q_acc += (size + dq) if qStrand == '+' else -(size + dq)
            return None  # 在该 chain 的 gap 内（被别的 chain 接住的情况由外层覆盖检查兜底）
        if maxend[j] <= pos0:
            return None  # 左侧所有 chain 均不覆盖
        j -= 1
    return None


def convert_bed(idx, in_bed, out_bed, unmapped):
    n_ok = n_bad = 0
    with open(in_bed) as fi, open(out_bed, 'w') as fo, open(unmapped, 'w') as fu:
        for line in fi:
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            p = line.split('\t')
            try:
                chrom = p[0]; start = int(p[1]); end = int(p[2])
            except (ValueError, IndexError):
                fu.write(line + '\tBAD_FORMAT\n'); n_bad += 1; continue
            r1 = convert_point(idx, chrom, start)
            r2 = convert_point(idx, chrom, max(start, end - 1))
            if r1 is None or r2 is None or r1[0] != r2[0]:
                fu.write(line + '\tUNMAPPED\n'); n_bad += 1; continue
            q1, p1 = r1; _, p2 = r2
            if p1 > p2: p1, p2 = p2, p1
            p2 += 1  # 恢复端点
            fo.write('\t'.join([q1, str(p1), str(p2)] + p[3:]) + '\n')
            n_ok += 1
    return n_ok, n_bad

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'build_index':
        idx = build_index(sys.argv[2])
        with open(sys.argv[3], 'wb') as f:
            pickle.dump(idx, f, protocol=4)
        print('chains indexed:', sum(len(v) for k, v in idx.items() if not k.endswith('__starts')))
    elif cmd == 'convert':
        with open(sys.argv[2], 'rb') as f:
            idx = pickle.load(f)
        ok, bad = convert_bed(idx, sys.argv[3], sys.argv[4], sys.argv[5])
        print(f'converted={ok} unmapped={bad}')
