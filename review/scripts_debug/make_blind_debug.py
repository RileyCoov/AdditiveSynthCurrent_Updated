"""Build review/blind_debug/ (neutral names) + review/.blind_key_debug.tsv from the Xcode Debug renders.
Same 17 comparisons as work2/fixes/scripts/make_blind.py. A fresh seed is drawn from os.urandom and written ONLY
to the key file; nothing about side assignment is printed."""
import random, shutil, os, hashlib, subprocess
R = '/Users/Riley/Downloads/AdditiveSynthCurrent-LongWindowFreqUsedForShortWindows/review'; B = f'{R}/blind_debug'; KEY = f'{R}/.blind_key_debug.tsv'
assert not os.path.exists(KEY), 'key already exists; refusing to overwrite'
os.makedirs(B, exist_ok=True)
PAIRS = [  # (file tag, version, reference cfg, test cfg) -- copied from the review's make_blind.py
 ('1985', 'dn5', 'base', 'cmb'), ('1985', 'dn5', 'base', 'f1'), ('drum', 'up5', 'base', 'f1'), ('drum', 'up5', 'base', 'cmb'),
 ('1985', 'up5', 'base', 'cmb'), ('fair', 'up5', 'base', 'f2'), ('fair', 'dn5', 'base', 'f2'), ('1985', 'dn5', 'f1', 'cmb'),
 ('1985', 'dn5', 'base', 'f3'), ('1985', 'up5', 'base', 'f3'), ('drum', 'dn5', 'base', 'f1'), ('happy', 'dn5', 'base', 'cmb'),
 ('fair', 'up5', 'base', 'f1'), ('1985', 'dno', 'base', 'cmb'), ('1985', 'dn5', 'base', 'f1b'), ('drum', 'dn5', 'base', 'f1b'),
 ('piano', 'dn5', 'base', 'cmb')]
CFG = {'base': 'baseline (all switches off)', 'f1': 'F1 SHORT_BASS_AMP_HOLD_HZ=400', 'f1b': 'F1b TSA_SPLIT_HZ=400', 'f2': 'F2 SHORT_BLOCK_8=1',
       'f3': 'F3 ATOM_MERGE_BINS=1.0', 'cmb': 'F1+F2 SHORT_BASS_AMP_HOLD_HZ=400 SHORT_BLOCK_8=1'}
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
seed = int.from_bytes(os.urandom(8), 'big'); rng = random.Random(seed)
commit = subprocess.run(['git', '-C', R, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
rows = [f'# seed\t{seed}', f'# binary\treview/bin_debug/patched md5 {md5(R + "/bin_debug/patched")}', f'# branch HEAD at pairing\t{commit}',
        'pair\tfile\tversion\tA\tB\tA_md5\tB_md5']
for i, (f, v, ref, tst) in enumerate(PAIRS, 1):
    a, b = (ref, tst) if rng.random() < 0.5 else (tst, ref); md = []
    for side, c in (('A', a), ('B', b)):
        src = f'{R}/renders_debug/{c}_{f}_{v}.wav'; dst = f'{B}/pair{i:02d}_{side}.wav'
        shutil.copyfile(src, dst); md.append(md5(dst)); assert md[-1] == md5(src)
    assert md[0] != md[1], f'pair{i:02d}: both sides identical'
    rows.append(f'pair{i:02d}\t{f}\t{v}\t{CFG[a]}\t{CFG[b]}\t{md[0]}\t{md[1]}')
open(KEY, 'w').write('\n'.join(rows) + '\n'); os.chmod(KEY, 0o600)
print(len(PAIRS), 'pairs written; every pair has two distinct renders; key written (not shown)')
