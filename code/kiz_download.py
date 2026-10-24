#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CattleGTEx mammary eQTL downloader (chunked Range resume). Usage: python kiz_download.py"""
import urllib.request, sys, ssl, os, time
import urllib.parse
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
URL = 'https://file.kiz.ac.cn/cattlegtex/' + urllib.parse.quote('Molecular QTLs/Single-Tissue cis-QTL Data/eQTL/Mammary_Gland.tar.gz')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'raw', 'cattlegtex_Mammary_Gland_eQTL.tar.gz')
TOTAL = 1339313499
CH = 1 << 23  # 8MB/段（服务器 69% 后频繁重置，缩小分块）
RESET_WAIT = 20  # 连续重置时的退避秒数
while True:
    have = os.path.getsize(OUT) if os.path.exists(OUT) else 0
    if have >= TOTAL:
        print('COMPLETE', have); break
    end = min(have + CH, TOTAL) - 1
    req = urllib.request.Request(URL, headers={'User-Agent': 'Mozilla/5.0', 'Range': f'bytes={have}-{end}'})
    try:
        with urllib.request.urlopen(req, timeout=120, context=ctx) as r, open(OUT, 'ab') as f:
            got = r.read()
            if not got:
                print('empty at', have); break
            f.write(got)
    except Exception as e:
        print('chunk err at', have, repr(e)[:80], flush=True); time.sleep(RESET_WAIT)
    print(f'  {os.path.getsize(OUT)/TOTAL*100:5.1f}%  {os.path.getsize(OUT)/1e6:.0f}MB', flush=True)
