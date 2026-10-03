from common import *
n='1985'
xo,fs=load(n,'original')
for name,(lo,hi) in {'kick(30-150)':(30,150),'snare(150-1500)':(150,1500),'hat(6k-16k)':(6000,16000),'broad(30-16k)':(30,16000)}.items():
    o,_=onsets(xo,fs,lo,hi)
    print(name, len(o), np.round(o,3).tolist())
