#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dryad ranged resume: the file_stream endpoint honours Range requests (206 verified); resumes from the local byte count.
Usage: python dryad_resume.py <file_stream_id> <local file path>
"""
import ssl, sys, os, re, json, time, http.cookiejar, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dryad_download as dd

def resume(fid, out_path, max_try=8):
    L = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    for attempt in range(max_try):
        try:
            req = urllib.request.Request(f'https://datadryad.org/downloads/file_stream/{fid}',
                                         headers={**dd.UA, 'Range': f'bytes={L}-'})
            r = dd.OP.open(req, timeout=180)
            first = r.read(512)
            if b'anubis' in first or b'Validating' in first:
                body = first + r.read(200000)
                dd.solve(body.decode('utf-8', errors='replace'), f'/downloads/file_stream/{fid}')
                time.sleep(1); continue
            if r.status != 206:
                print('server ignored Range, status', r.status, file=sys.stderr)
                return None
            with open(out_path, 'ab') as f:
                f.write(first); tot = len(first)
                while True:
                    c = r.read(1 << 20)
                    if not c: break
                    f.write(c); tot += len(c)
            L += tot
            print(f'attempt {attempt}: +{tot}B, total {L}', file=sys.stderr)
            # 检查是否完成：再发一次 Range 请求末尾 1 字节
            req2 = urllib.request.Request(f'https://datadryad.org/downloads/file_stream/{fid}',
                                          headers={**dd.UA, 'Range': f'bytes={L-1}-{L}'})
            r2 = dd.OP.open(req2, timeout=60)
            r2.read()
            cr = r2.headers.get('Content-Range', '')
            m = re.search(r'/(\d+)$', cr)
            if m and int(m.group(1)) <= L:
                print('COMPLETE', L, file=sys.stderr)
                return L
            time.sleep(2)
        except Exception as e:
            print('retry:', repr(e)[:120], file=sys.stderr); time.sleep(4)
            L = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    return None

if __name__ == '__main__':
    n = resume(sys.argv[1], sys.argv[2])
    print('OK' if n else 'INCOMPLETE', n)
