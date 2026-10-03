import pickle, numpy as np, csv
R={'up_fifth':2**(7/12),'down_fifth':2**(-7/12)}
ev=pickle.load(open('an/1985_events.pkl','rb'))
rows=[]
for v in ['up_fifth','down_fifth']:
    other='down_fifth' if v=='up_fifth' else 'up_fifth'
    for e in ev[v]:
        if e['rel_to_frame']<=-12: continue
        fs_=e['f']/R[v]
        twin=any(abs(fs_-b['f']/R[other])/fs_<0.03 and e['t0']<=b['t1']+0.02 and b['t0']<=e['t1']+0.02 for b in ev[other])
        rows.append([v,round(e['t0'],3),round(e['t1'],3),int(e['t0']*48000),int(e['t1']*48000),round(e['f'],1),round(fs_,1),e['nfr'],round(e['exc'],1),round(e['rel_to_frame'],1),twin])
rows.sort(key=lambda r:(r[1],r[0]))
with open('deliver/1985_prominent_tonal_events.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['output','t0_s(window centre)','t1_s','sample0','sample1','output_freq_Hz','source_freq_Hz','hops(256)','excess_over_mapped_source_dB','level_re_frame_max_dB','twin_in_other_direction'])
    w.writerows(rows)
print(len(rows),'rows; twins:',sum(r[-1] for r in rows))
