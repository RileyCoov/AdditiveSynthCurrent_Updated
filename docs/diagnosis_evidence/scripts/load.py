import numpy as np, soundfile as sf, os
# DIAG_ROOT: folder holding original/ unity/ up_fifth/ down_fifth/ (default = the shipped engine state)
ROOT=os.environ.get('DIAG_ROOT', os.path.expanduser('~/Desktop/comparingOut/ENGINE_STATE_2026-10-02'))
VERS=['original','unity','up_fifth','down_fifth']
NAMES=sorted([f[:-len('_original.wav')] for f in os.listdir(f'{ROOT}/original') if f.endswith('.wav')])
def load(name, ver):
    x, fs = sf.read(f'{ROOT}/{ver}/{name}_{ver}.wav', always_2d=True)
    return x.mean(1).astype(np.float64), fs
def info(name, ver):
    return sf.info(f'{ROOT}/{ver}/{name}_{ver}.wav')
