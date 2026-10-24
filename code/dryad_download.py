#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dryad downloader with automatic Anubis proof-of-work handling.
Since 2026-09 Dryad serves an Anubis challenge on /downloads/file_stream/* and parts of the API:
  the challenge page embeds JSON {rules:{difficulty:N}, challenge:{id, issuedAt, randomData, ...}}
  protocol (fast algorithm): find a nonce such that sha256(randomData + str(nonce)) starts with N '0's
  pass: GET /.within.website/x/cmd/anubis/api/pass-challenge?id=&response=&nonce=&redir=&elapsedTime=200
  once the techaro.lol-anubis-auth cookie is set, the original direct link streams the file.
Usage:
  python dryad_download.py <dataset_doi_tail>            # list file_stream ids of the dataset
  python dryad_download.py <dataset_doi_tail> <file_id> <output_path>
"""
import ssl, sys, http.cookiejar, urllib.request, re, json, hashlib, time
from urllib.parse import urlencode

CTX = ssl.create_default_context(); CTX.check_hostname=False; CTX.verify_mode=ssl.CERT_NONE
UA = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36','Accept':'*/*'}
BASE = 'https://datadryad.org'
_cj = http.cookiejar.CookieJar()
OP = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(_cj), urllib.request.HTTPSHandler(context=CTX))

def solve(body, redir):
    m = re.search(r'<script[^>]*>\s*({.*?"difficulty".*?})\s*</script>', body, re.S)
    if not m: return False
    cfg = json.loads(m.group(1)); ch = cfg['challenge']; diff = cfg['rules']['difficulty']
    base = ch['randomData']
    n, h = None, None
    for i in range(50_000_000):
        hh = hashlib.sha256(f'{base}{i}'.encode()).hexdigest()
        if hh[:diff] == '0'*diff: n, h = i, hh; break
    if n is None: return False
    q = urlencode({'id':ch['id'], 'response':h, 'nonce':n, 'redir':redir, 'elapsedTime':'200'})
    OP.open(urllib.request.Request(f'{BASE}/.within.website/x/cmd/anubis/api/pass-challenge?{q}', headers=UA), timeout=60).read()
    return True

def stream(fid, out_path, max_try=4):
    for _ in range(max_try):
        try:
            r = OP.open(urllib.request.Request(f'{BASE}/downloads/file_stream/{fid}', headers=UA), timeout=120)
            first = r.read(512)
            if b'anubis' in first or b'Validating' in first:
                body = first + r.read(200000)
                solve(body.decode('utf-8', errors='replace'), f'/downloads/file_stream/{fid}')
                time.sleep(1); continue
            tot = 0
            with open(out_path, 'wb') as f:
                f.write(first); tot = len(first)
                while True:
                    c = r.read(1<<20)
                    if not c: break
                    f.write(c); tot += len(c)
            return tot
        except Exception as e:
            print('retry:', repr(e)[:120], file=sys.stderr); time.sleep(3)
    return None

def list_ids(doi_tail):
    body = OP.open(urllib.request.Request(f'{BASE}/dataset/doi:{doi_tail}', headers=UA), timeout=60).read().decode('utf-8', errors='replace')
    return list(dict.fromkeys(re.findall(r'/downloads/file_stream/(\d+)', body)))

if __name__ == '__main__':
    if len(sys.argv) == 2:
        print('\n'.join(list_ids(sys.argv[1])))
    elif len(sys.argv) == 4:
        n = stream(sys.argv[2], sys.argv[3])
        print('OK' if n is not None else 'FAIL', n)
    else:
        print(__doc__)
