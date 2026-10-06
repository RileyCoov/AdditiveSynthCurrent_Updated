"""Pre-registered metrics for the fix prototypes (definitions frozen in fix_prototypes.md, 2026-10-05T13:53:58-07:00).
usage: fixmetrics.py <tag> [<tag> ...]   (tag = <cfg>_<file>_<ver>, render renders/<tag>.wav, log logs/<tag>.log)
       fixmetrics.py --table
M1 reuses Analyst B's sig.py unchanged; HFATT reuses Analyst A's s09 definition; CLICK reuses E's s19 definition;
BIRDIES reuses A's s16 method. OV and RISE are defined in the pre-registration."""
import sys, json, numpy as np, soundfile as sf
from scipy.signal import butter, sosfiltfilt, fftconvolve
from scipy.ndimage import label, maximum_filter1d
sys.path.insert(0, '/Users/Riley/Downloads/work2/1985_causal/scripts'); import sig            # B
sys.path.insert(0, '/Users/Riley/Downloads/work2/1985_signal'); import common as ac           # A
ac.WORK = '/Users/Riley/Downloads/work2/1985_signal'   # PATH ONLY (Debug run): RB3 control location for load_rb
F = '/Users/Riley/Downloads/AdditiveSynthCurrent-LongWindowFreqUsedForShortWindows/review'; A = '/Users/Riley/Desktop/comparingOut/EXPERT_REVIEW_2026-10-04/audio/original'; SR = 48000
NAMES = {'1985': '1985', 'drum': 'DrumLoopShort', 'fair': 'FairlightAdditiveWaveTA1_C2', 'piano': 'PianoSampleMono', 'happy': 'HappyMono'}
SEMI = {'dn5': -7, 'up5': 7, 'dno': -12, 'uni': 0}
RBV = {'dn5': 'down_fifth', 'up5': 'up_fifth'}

def ov(x, y, r, Ss, lo, hi):
    idx, ex = sig.env(sig.band(x, lo, hi)); _, ey = sig.env(sig.band(y, lo * r, hi * r)); D = 20 * np.log10(ey / ex); v = []
    for S in Ss:
        pre = (idx >= S - 2048) & (idx < S + 512); sb = (idx >= S + 1472) & (idx < S + 2496)
        if pre.sum() < 10 or sb.sum() < 10: continue
        v.append(np.mean(D[sb] - np.median(D[pre])))
    return float(np.median(v)) if v else None

def env5(x):
    w = np.hanning(241); w /= w.sum(); return np.sqrt(np.maximum(fftconvolve(x * x, w, 'same'), 1e-20))[::8]

def t50(w):
    k = int(np.argmax(w)); P = w[k]; j = k
    while j > 0 and w[j - 1] >= 0.5 * P: j -= 1
    return j

def rise(x, y, r, Ss, lo, hi):
    ex = env5(sig.band(x, lo, hi)); ey = env5(sig.band(y, lo * r, hi * r)); out = []
    for S in Ss:
        a = S // 8; b = min(len(ex), (S + 6144) // 8)
        if b - a < 200: continue
        # pre-level window [S-2048,S); for S<2048 (block at file start) [0,S+512) -- deviation noted in the outcomes
        pre = np.median(ex[a - 256:a]) if a >= 256 else np.median(ex[0:a + 64])
        wx = ex[a:b]
        if 20 * np.log10(wx.max() / pre) < 10: continue
        jo = t50(wx); jy = t50(ey[a:b]); out.append(dict(S=int(S), orig=jo * 8, d_ms=(jy - jo) * 8 / 48.0))
    return out

def hfatt(x, Ss):            # A's s09, unchanged
    y = ac.bandpass_fft(x, 6000, 16000); e = 20 * np.log10(np.sqrt(np.convolve(y ** 2, np.ones(24) / 24, 'same'))[::4] + 1e-7)
    res = []
    for S in Ss:
        a = S // 4; w = e[a + 512 // 4:a + 3584 // 4]
        if len(w) < 700: continue
        sl = np.convolve(np.diff(w), np.ones(60) / 60, 'same'); res.append(512 + int(np.argmax(sl)) * 4)
    return res

def click(x, Ss):            # E's s19, unchanged
    x = x * 32768; d2 = np.abs(np.diff(x, 2)); hp = sosfiltfilt(butter(4, 4000, btype='highpass', fs=SR, output='sos'), x); a = []; c = []
    for S in Ss:
        if S + 2600 > len(x): continue
        a.append(d2[S + 2488:S + 2504].max())
        c.append(20 * np.log10(np.abs(hp[S + 2480:S + 2512]).max() / (np.sqrt(np.mean(hp[S + 2048:S + 2400] ** 2)) + 1e-9)))
    return float(np.median(c)), float(np.median(a)), len(c)

N = 2048; H = 256; WN = np.hanning(N); FQ = np.fft.rfftfreq(N, 1 / SR)
def stftdb(x):
    nf = (len(x) - N) // H + 1; out = np.empty((nf, N // 2 + 1), np.float32)
    for i in range(nf): out[i] = 20 * np.log10(np.abs(np.fft.rfft(x[i * H:i * H + N] * WN)) + 1e-9)
    return out
def blobs(De, Dref, Do, r, band):  # A's s16, unchanged logic (count only)
    n = min(len(De), len(Dref), len(Do)); De, Dref, Do = De[:n], Dref[:n], Do[:n]
    Dom = np.empty_like(De)
    for i in range(n): Dom[i] = np.interp(FQ / r, FQ, Do[i], right=-200)
    E = De - np.maximum(Dref, Dom); fmax = De.max(axis=1, keepdims=True); fs = FQ / r; bm = (fs >= band[0]) & (fs < band[1])
    M = (E >= 10) & (De >= fmax - 50) & bm[None, :]; lab, k = label(M, structure=np.ones((3, 3)))
    sizes = np.bincount(lab.ravel())[1:]; return int((sizes >= 2).sum())

def run(tag):
    cfg, f, ver = tag.split('_'); name = NAMES[f]; semi = SEMI[ver]; r = 2 ** (semi / 12)
    x = sig.load(f'{A}/{name}_original.wav'); y = sig.load(f'{F}/renders_debug/{tag}.wav'); n = min(len(x), len(y)); x, y = x[:n], y[:n]
    Ss = sig.segments_from_log(f'{F}/logs_debug/{tag}.log'); o = dict(tag=tag, cfg=cfg, file=f, ver=ver, nseg=len(Ss))
    if semi == 0: return o
    o['lag'] = sig.lag_check(x, y, r)
    res, _, _ = sig.m1(x, y, r, Ss)
    if res:
        o['M1_dip'] = float(np.median([q['dip'] for q in res])); o['M1_rms'] = float(np.median([q['rms'] for q in res])); o['M1_n'] = len(res)
        for lo, hi in [(30, 160), (160, 400), (400, 2000)]: o[f'OV_{lo}_{hi}'] = ov(x, y, r, Ss, lo, hi)
    if f in ('1985', 'drum', 'happy'):
        rr = rise(x, y, r, Ss, 40, 160); o['RISE_kick_ms'] = float(np.median([q['d_ms'] for q in rr])) if rr else None; o['RISE_kick_n'] = len(rr)
        o['RISE_kick_per'] = [round(q['d_ms'], 1) for q in rr]
    if f == 'piano':
        for lo, hi, k in [(80, 400, 'LF'), (2000, 8000, 'HF')]:
            rr = rise(x, y, r, Ss, lo, hi); o[f'RISE_{k}_ms'] = [round(q['d_ms'], 2) for q in rr]; o[f'RISE_{k}_S'] = [q['S'] for q in rr]
    if f == '1985':
        h = hfatt(y, Ss); ho = hfatt(x, Ss); o['HFATT_med'] = float(np.median(h)); o['HFATT_orig_med'] = float(np.median(ho))
        o['HFATT_absdev_med'] = float(np.median(np.abs(np.array(h) - np.array(ho))))
        if ver in RBV:
            De = stftdb(y); Dr = stftdb(ac.load_rb(RBV[ver])); Do = maximum_filter1d(stftdb(x), 5, axis=1); dur = n / SR
            o['BIRD_400_3500_per_s'] = blobs(De, Dr, Do, r, (400, 3500)) / dur; o['BIRD_40_400'] = blobs(De, Dr, Do, r, (40, 400))
    if f in ('fair', 'drum'):
        o['CLICK_db'], o['CLICK_d2'], o['CLICK_n'] = click(y, Ss); o['CLICK_orig_db'] = click(x, Ss)[0]
    return o

def table():
    rows = [json.loads(l) for l in open(f'{F}/logs_debug/metrics.jsonl')]; last = {}
    for q in rows: last[q['tag']] = q
    keys = ['lag', 'M1_dip', 'M1_rms', 'OV_30_160', 'OV_160_400', 'OV_400_2000', 'RISE_kick_ms', 'RISE_kick_n', 'HFATT_med', 'HFATT_absdev_med',
            'BIRD_400_3500_per_s', 'BIRD_40_400', 'CLICK_db', 'CLICK_d2', 'CLICK_orig_db', 'RISE_LF_ms', 'RISE_HF_ms']
    for fv in sorted({(q['file'], q['ver']) for q in last.values() if q['ver'] != 'uni'}):
        print(f'== {fv[0]} {fv[1]}')
        for q in sorted([q for q in last.values() if (q['file'], q['ver']) == fv], key=lambda q: (q['cfg'] != 'base', q['cfg'])):
            print('  %-6s' % q['cfg'] + ' '.join(f'{k}={q[k]:.2f}' if isinstance(q.get(k), float) else f'{k}={q[k]}' for k in keys if k in q and q[k] is not None))

if __name__ == '__main__':
    if sys.argv[1] == '--table': table(); sys.exit()
    with open(f'{F}/logs_debug/metrics.jsonl', 'a') as fh:
        for t in sys.argv[1:]:
            o = run(t); fh.write(json.dumps(o) + '\n'); fh.flush(); print(json.dumps(o))
