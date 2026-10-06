"""Score Debug metrics (logs_debug/metrics.jsonl) against the frozen predictions in predictions_debug.md.
Release values and tolerances are copied from that file's frozen section. Prints the OUTCOMES markdown."""
import json
F = '/Users/Riley/Downloads/AdditiveSynthCurrent-LongWindowFreqUsedForShortWindows/review'
M = {}
for l in open(f'{F}/logs_debug/metrics.jsonl'):
    q = json.loads(l); M[q['tag']] = q
g = lambda c, f, v, k: M[f'{c}_{f}_{v}'][k]
V3 = ['dn5', 'up5', 'dno']
rows, flags = [], []
def fmt(x, nd=2):
    return '/'.join(f'{v:+.{nd}f}' for v in x) if isinstance(x, (list, tuple)) else f'{x:+.{nd}f}'
def cmp_rel(label, deb, rel, tol, unit='', pct=False, minabs=0.0):
    for d, r in zip(deb if isinstance(deb, list) else [deb], rel if isinstance(rel, list) else [rel]):
        t = max(abs(r) * tol, minabs) if pct else tol
        if abs(d - r) > t + 1e-9:
            flags.append(f'| {label} | {r:g} | {d:g} | {d - r:+.2f}{unit} | ±{t:g}{unit} |')
def item(fix, name, pred, deb_txt, interval, rule):
    rows.append(f'| {fix} | {name} | {pred} | {deb_txt} | {interval} | {rule} |')
inr = lambda x, lo, hi: lo - 1e-9 <= x <= hi + 1e-9
hitall = lambda xs, ps, tols: all(abs(x - p) <= t + 1e-9 for x, p, t in zip(xs, ps, tols))

# ---- F1
b = [g('base', '1985', v, 'M1_dip') for v in V3]; x = [g('f1', '1985', v, 'M1_dip') for v in V3]
cmp_rel('base 1985 M1_dip dn5/up5/dno', b, [-5.61, -4.96, -4.81], 0.6, ' dB'); cmp_rel('F1 1985 M1_dip', x, [-1.40, -1.19, -0.56], 0.6, ' dB')
item('F1', 'P1 1985 M1_dip dn5/up5/dno', '−1.4/−1.2/−0.6 ±0.6', f'{fmt(x)} (base {fmt(b)})',
     'hit' if hitall(x, [-1.4, -1.2, -0.6], [.6] * 3) else 'miss', 'PASS' if all(inr(v, -2, 1) for v in x) else '**FAIL**')
b = [g('base', '1985', v, 'OV_160_400') for v in V3]; x = [g('f1', '1985', v, 'OV_160_400') for v in V3]
cmp_rel('base 1985 OV160-400', b, [-4.25, -4.50, -3.76], 1.0, ' dB'); cmp_rel('F1 1985 OV160-400', x, [-2.47, -1.25, 0.12], 1.0, ' dB')
miss = [v for v, xx in zip(V3, x) if abs(xx + 1.0) > 1.0 + 1e-9]
item('F1', 'P2 1985 OV 160-400 dn5/up5/dno', '−1.0 ±1.0', f'{fmt(x)} (base {fmt(b)})', 'hit' if not miss else 'miss ' + ', '.join(miss),
     'PASS' if all(inr(v, -2.5, 1.5) for v in x) else '**FAIL**')
b = [g('base', '1985', v, 'OV_400_2000') for v in V3]; x = [g('f1', '1985', v, 'OV_400_2000') for v in V3]; d = [xx - bb for xx, bb in zip(x, b)]
cmp_rel('base 1985 OV400-2000', b, [-4.82, -4.56, -4.08], 0.7, ' dB'); cmp_rel('F1 1985 OV400-2000 Δ vs base', d, [0.35, 0.11, 0.09], 0.7, ' dB')
item('F1', 'P2 1985 OV 400-2000 Δ vs base', '±0.7', fmt(d), 'hit' if all(abs(v) <= .7 for v in d) else 'miss',
     'PASS' if all(v >= -1.0 for v in d) else '**FAIL**')
bu, bd = g('base', 'drum', 'up5', 'RISE_kick_ms'), g('base', 'drum', 'dn5', 'RISE_kick_ms')
du, dd = g('f1', 'drum', 'up5', 'RISE_kick_ms') - bu, g('f1', 'drum', 'dn5', 'RISE_kick_ms') - bd
cmp_rel('base Drum kick RISE up5/dn5', [bu, bd], [11.33, 32.67], 3.0, ' ms'); cmp_rel('F1 Drum kick RISE Δ up5/dn5', [du, dd], [1.33, 0.0], 3.0, ' ms')
item('F1', 'P3 Drum kick RISE Δ up5/dn5', '+3 ±3 ms; "expect FAIL"', f'{fmt([du, dd])} ms (base {bu:+.2f}/{bd:+.2f})',
     'hit' if all(abs(v - 3) <= 3 for v in (du, dd)) else 'miss', 'PASS' if du <= 2 and dd <= 2 else '**FAIL**')
b, x = g('base', 'happy', 'dn5', 'M1_dip'), g('f1', 'happy', 'dn5', 'M1_dip')
cmp_rel('base HappyMono dn5 M1_dip', b, -4.75, 0.5, ' dB'); cmp_rel('F1 HappyMono dn5 M1_dip', x, -0.85, 0.5, ' dB')
item('F1', 'P4 HappyMono dn5 M1_dip', '[−2.5, +1.0]', f'{x:+.2f} (base {b:+.2f})', 'hit' if inr(x, -2.5, 1) else 'miss', 'PASS' if inr(x, -2.5, 1) else '**FAIL**')
bd_, bu_ = g('base', 'fair', 'dn5', 'CLICK_db'), g('base', 'fair', 'up5', 'CLICK_db'); xd, xu = g('f1', 'fair', 'dn5', 'CLICK_db'), g('f1', 'fair', 'up5', 'CLICK_db')
cmp_rel('base FairC2 CLICK dn5/up5', [bd_, bu_], [17.27, 20.13], 3.0, ' dB'); cmp_rel('F1 FairC2 CLICK dn5/up5', [xd, xu], [10.14, 15.94], 3.0, ' dB')
ok = abs(xd - bd_) <= 3 and abs(xu - bu_) <= 3
item('F1', 'P5 FairC2 CLICK within ±3 dB of base', '±3 dB', f'dn5 {bd_:.2f} → {xd:.2f}; up5 {bu_:.2f} → {xu:.2f}', 'hit' if ok else 'miss (lower)', 'PASS' if ok else '**FAIL**')
pb = (g('base', 'piano', 'dn5', 'RISE_LF_ms')[0], g('base', 'piano', 'dn5', 'RISE_HF_ms')[0])
p1 = (g('f1', 'piano', 'dn5', 'RISE_LF_ms')[0] - pb[0], g('f1', 'piano', 'dn5', 'RISE_HF_ms')[0] - pb[1])
cmp_rel('F1 Piano RISE Δ LF/HF', list(p1), [0.0, 0.0], 3.0, ' ms')
item('F1', 'P6 Piano dn5 RISE Δ LF/HF', '±3 ms', f'{fmt(list(p1))} ms', 'hit' if all(abs(v) <= 3 for v in p1) else 'miss', 'PASS' if all(v <= 2 for v in p1) else '**FAIL**')
d = [g('f1', '1985', v, 'HFATT_med') - g('base', '1985', v, 'HFATT_med') for v in V3]
cmp_rel('base 1985 HFATT dn5/up5', [g('base', '1985', v, 'HFATT_med') for v in V3[:2]], [2284, 2440], 150, ' smp')
cmp_rel('F1 1985 HFATT Δ', d, [0, 0, 0], 150, ' smp')
item('F1', 'P7 1985 HFATT Δ dn5/up5/dno', '±150 smp', '/'.join(f'{v:+.0f}' for v in d), 'hit' if all(abs(v) <= 150 for v in d) else 'miss', 'PASS' if all(abs(v) <= 150 for v in d) else '**FAIL**')
d = [g('f1', '1985', v, 'RISE_kick_ms') - g('base', '1985', v, 'RISE_kick_ms') for v in V3]
cmp_rel('F1 1985 kick RISE Δ', d, [-3.4, 0.7, 1.3], 2.0, ' ms')
bad = [v for v, dd_ in zip(V3, d) if dd_ < -1 - 1e-9]
item('F1', 'P7 1985 kick RISE Δ not earlier than −1 ms', '—', fmt(d, 1) + ' ms', '—', 'PASS' if not bad else f'**FAIL ({", ".join(bad)})**')

# ---- F1+F2
x = [g('cmb', '1985', v, 'M1_dip') for v in V3]; f1 = [g('f1', '1985', v, 'M1_dip') for v in V3]
cmp_rel('F1+F2 1985 M1_dip', x, [-0.89, -0.51, -1.94], 0.7, ' dB')
miss = [v for v, a, c in zip(V3, x, f1) if abs(a - c) > .7 + 1e-9]
item('F1+F2', 'M1_dip vs F1 ±0.7; rule [−2, +1]', '—', fmt(x), 'hit' if not miss else 'miss ' + ', '.join(miss), 'PASS' if all(inr(v, -2, 1) for v in x) else '**FAIL**')
x = [g('cmb', '1985', v, 'OV_160_400') for v in V3]; f1 = [g('f1', '1985', v, 'OV_160_400') for v in V3]
cmp_rel('F1+F2 1985 OV160-400', x, [-3.43, -1.26, -2.66], 0.7, ' dB')
miss = [v for v, a, c in zip(V3, x, f1) if abs(a - c) > .7 + 1e-9]; bad = [v for v, a in zip(V3, x) if not inr(a, -2.5, 1.5)]
item('F1+F2', 'OV 160-400 vs F1 ±0.7; rule [−2.5, +1.5]', '—', fmt(x), 'hit' if not miss else 'miss ' + ', '.join(miss), 'PASS' if not bad else f'**FAIL ({", ".join(bad)})**')
cf = [g('cmb', 'fair', 'up5', 'CLICK_db'), g('cmb', 'fair', 'dn5', 'CLICK_db')]; cdr = [g('cmb', 'drum', 'up5', 'CLICK_db'), g('cmb', 'drum', 'dn5', 'CLICK_db')]
f2f = [g('f2', 'fair', 'up5', 'CLICK_db'), g('f2', 'fair', 'dn5', 'CLICK_db')]; f2d = [g('f2', 'drum', 'up5', 'CLICK_db'), g('f2', 'drum', 'dn5', 'CLICK_db')]
bf = [g('base', 'fair', 'up5', 'CLICK_db'), g('base', 'fair', 'dn5', 'CLICK_db')]
cmp_rel('F1+F2 CLICK Fair up5/dn5', cf, [5.61, 3.52], 2.0, ' dB'); cmp_rel('F1+F2 CLICK Drum up5/dn5', cdr, [9.45, 7.82], 2.0, ' dB')
hit = all(abs(a - c) <= 2 for a, c in zip(cf + cdr, f2f + f2d)); rule = all(b_ - a >= 6 and a <= 10 for a, b_ in zip(cf, bf))
item('F1+F2', 'CLICK vs F2 ±2 dB', '—', f'Fair up5 {cf[0]:.2f}, dn5 {cf[1]:.2f}; Drum {cdr[0]:.2f}/{cdr[1]:.2f}', 'hit' if hit else 'miss', 'PASS' if rule else '**FAIL**')
ku, kd = g('cmb', 'drum', 'up5', 'RISE_kick_ms'), g('cmb', 'drum', 'dn5', 'RISE_kick_ms')
cmp_rel('F1+F2 Drum kick RISE up5/dn5', [ku, kd], [25.00, 27.17], 1.0, ' ms')
fu, fd = g('f1', 'drum', 'up5', 'RISE_kick_ms'), g('f1', 'drum', 'dn5', 'RISE_kick_ms')
bad = [n for n, a, b_ in (('up5', ku, bu), ('dn5', kd, bd)) if a > b_ + 2 + 1e-9]
item('F1+F2', 'Drum kick RISE vs F1 ±1 ms; rule ≤ base + 2 ms', '—', f'up5 {ku:.2f}, dn5 {kd:.2f} ms (base {bu:+.2f}/{bd:+.2f})',
     'hit' if abs(ku - fu) <= 1 and abs(kd - fd) <= 1 else 'miss', 'PASS' if not bad else f'**FAIL ({", ".join(bad)})**')

# ---- F2
xu, xd = g('f2', 'fair', 'up5', 'CLICK_db'), g('f2', 'fair', 'dn5', 'CLICK_db')
cmp_rel('F2 FairC2 CLICK up5/dn5', [xu, xd], [5.79, 3.48], 2.0, ' dB')
cmp_rel('F2 FairC2 CLICK max|d2| up5/dn5', [g('f2', 'fair', 'up5', 'CLICK_d2'), g('f2', 'fair', 'dn5', 'CLICK_d2')], [270, 76], 0.25, ' LSB', pct=True)
cmp_rel('base FairC2 CLICK max|d2| up5/dn5', [g('base', 'fair', 'up5', 'CLICK_d2'), g('base', 'fair', 'dn5', 'CLICK_d2')], [1403, 959], 0.25, ' LSB', pct=True)
ok = bu_ - xu >= 6 and bd_ - xd >= 6 and xu <= 10 and xd <= 10
item('F2', 'P1 FairC2 CLICK up5/dn5', 'drop ≥6 dB, ≤ +10 dB', f'{bu_:.2f} → {xu:.2f}; {bd_:.2f} → {xd:.2f} (orig {g("base","fair","up5","CLICK_orig_db"):.2f}); max|d2| {g("base","fair","up5","CLICK_d2"):.0f} → {g("f2","fair","up5","CLICK_d2"):.0f}, {g("base","fair","dn5","CLICK_d2"):.0f} → {g("f2","fair","dn5","CLICK_d2"):.0f} LSB',
     'hit' if ok else 'miss', 'PASS' if ok else '**FAIL**')
x = [g('f2', '1985', v, 'M1_dip') for v in V3]; b = [g('base', '1985', v, 'M1_dip') for v in V3]
cmp_rel('F2 1985 M1_dip dn5/up5', x[:2], [-5.22, -5.97], 0.5, ' dB'); cmp_rel('F2 1985 M1_dip dno', x[2], -5.55, 0.7, ' dB')
item('F2', 'P2 1985 M1_dip dn5/up5/dno', '−5.3/−5.9 ±0.5; −5.0 ±0.7', fmt(x), 'hit' if hitall(x, [-5.3, -5.9, -5.0], [.5, .5, .7]) else 'miss',
     'PASS (no harm)' if all(a >= c - 1.2 for a, c in zip(x, b)) else '**FAIL**')
du2, dd2 = g('f2', 'drum', 'up5', 'RISE_kick_ms') - bu, g('f2', 'drum', 'dn5', 'RISE_kick_ms') - bd
cmp_rel('F2 Drum kick RISE Δ up5/dn5', [du2, dd2], [-8.8, -5.5], 2.0, ' ms')
item('F2', 'P3 Drum kick RISE Δ up5/dn5', '±2 ms', f'{fmt([du2, dd2], 1)} ms', 'hit' if abs(du2) <= 2 and abs(dd2) <= 2 else 'miss (earlier)', 'PASS' if du2 <= 2 and dd2 <= 2 else '**FAIL**')
d = [g('f2', '1985', v, 'HFATT_med') - g('base', '1985', v, 'HFATT_med') for v in V3]
cmp_rel('F2 1985 HFATT Δ', d, [-20, 16, 32], 150, ' smp')
item('F2', 'P4 1985 HFATT Δ dn5/up5/dno', '±150 smp', '/'.join(f'{v:+.0f}' for v in d), 'hit' if all(abs(v) <= 150 for v in d) else 'miss', 'PASS' if all(abs(v) <= 150 for v in d) else '**FAIL**')
p2 = (g('f2', 'piano', 'dn5', 'RISE_LF_ms')[0] - pb[0], g('f2', 'piano', 'dn5', 'RISE_HF_ms')[0] - pb[1])
cmp_rel('F2 Piano RISE Δ LF/HF', list(p2), [0.0, 3.0], 2.0, ' ms')
bad = [n for n, v in zip(('LF', 'HF'), p2) if abs(v) > 2 + 1e-9]
item('F2', 'P4 Piano RISE Δ LF/HF (n=1)', '±2 ms', f'{fmt(list(p2))} ms', 'hit' if not bad else 'miss ' + ', '.join(bad), 'PASS' if not bad else f'**FAIL ({", ".join(bad)})**')

# ---- F3
bb = [g('base', '1985', v, 'BIRD_400_3500_per_s') for v in ('dn5', 'up5')]; xb = [g('f3', '1985', v, 'BIRD_400_3500_per_s') for v in ('dn5', 'up5')]
cmp_rel('base 1985 BIRDIES dn5/up5', bb, [21.09, 10.00], 0.15, ' /s', pct=True); cmp_rel('F3 1985 BIRDIES dn5/up5', xb, [9.35, 3.91], 0.15, ' /s', pct=True)
rel = [100 * (x_ / b_ - 1) for x_, b_ in zip(xb, bb)]
item('F3', 'BIRDIES 400-3500 Hz dn5/up5', '−15 ±15%', f'{bb[0]:.2f} → {xb[0]:.2f} /s ({rel[0]:+.0f}%); {bb[1]:.2f} → {xb[1]:.2f} /s ({rel[1]:+.0f}%)',
     'hit' if all(-30 <= r <= 0 for r in rel) else 'miss (bigger)', 'PASS' if all(r <= -10 for r in rel) else '**FAIL**')
b4 = [g('base', '1985', v, 'BIRD_40_400') for v in ('dn5', 'up5')]; x4 = [g('f3', '1985', v, 'BIRD_40_400') for v in ('dn5', 'up5')]
cmp_rel('base 1985 40-400 blobs', b4, [4, 6], 0.2, '', pct=True, minabs=1); cmp_rel('F3 1985 40-400 blobs', x4, [5, 2], 0.2, '', pct=True, minabs=1)
bad = [v for v, a, c in zip(('dn5', 'up5'), x4, b4) if a > 1.2 * c]
item('F3', '40-400 Hz blob count not up >20%', '—', f'dn5 {b4[0]} → {x4[0]}; up5 {b4[1]} → {x4[1]}', '—', 'PASS' if not bad else f'**FAIL ({", ".join(bad)})**')
x = [g('f3', '1985', v, 'M1_dip') for v in ('dn5', 'up5')]; b = [g('base', '1985', v, 'M1_dip') for v in ('dn5', 'up5')]
cmp_rel('F3 1985 M1_dip dn5/up5', x, [-5.61, -5.04], 0.5, ' dB')
ok = all(abs(a - c) <= .5 for a, c in zip(x, b))
item('F3', 'M1_dip within ±0.5 dB of base', '—', f'{fmt(x)} (base {fmt(b)})', 'hit' if ok else 'miss', 'PASS' if ok else '**FAIL**')

# ---- review "other measured effects" (not pass/fail items), Release comparison only
for tag, k, rel_, tol in [('f1_1985_dn5', 'OV_30_160', -1.41, .7), ('f1_1985_up5', 'OV_30_160', .23, .7), ('f1_1985_dno', 'OV_30_160', -.10, .7),
                          ('base_1985_dn5', 'OV_30_160', -6.34, .7), ('base_1985_up5', 'OV_30_160', -5.08, .7), ('base_1985_dno', 'OV_30_160', -5.03, .7),
                          ('base_fair_dn5', 'OV_160_400', .15, .7), ('base_fair_up5', 'OV_160_400', -2.96, .7),
                          ('f1_fair_dn5', 'OV_160_400', 4.62, .7), ('f1_fair_up5', 'OV_160_400', 4.41, .7),
                          ('base_happy_dn5', 'OV_30_160', -3.85, .7), ('f1_happy_dn5', 'OV_30_160', 2.75, .7),
                          ('cmb_happy_dn5', 'M1_dip', -.88, .5), ('cmb_happy_dn5', 'OV_30_160', -.36, .7),
                          ('f3_1985_dn5', 'OV_400_2000', -3.96, .7), ('f3_1985_dn5', 'HFATT_absdev_med', 504, 150), ('base_1985_dn5', 'HFATT_absdev_med', 612, 150),
                          ('base_drum_dn5', 'CLICK_db', 13.98, 2), ('base_drum_up5', 'CLICK_db', 10.60, 2), ('f2_drum_dn5', 'CLICK_db', 7.83, 2), ('f2_drum_up5', 'CLICK_db', 9.45, 2),
                          ('f1b_1985_dn5', 'M1_dip', .64, 1.0), ('f1b_1985_dn5', 'OV_160_400', 2.40, 1.5)]:
    cmp_rel(f'{tag} {k} (other effect)', M[tag][k], rel_, tol)
for tag, k, rel_ in [('f1_fair_dn5', 'CLICK_d2', 393), ('f1_fair_up5', 'CLICK_d2', 821), ('base_drum_dn5', 'CLICK_d2', 968), ('base_drum_up5', 'CLICK_d2', 1699),
                     ('f2_drum_dn5', 'CLICK_d2', 154), ('f2_drum_up5', 'CLICK_d2', 726)]:
    cmp_rel(f'{tag} {k} (other effect)', M[tag][k], rel_, .25, ' LSB', pct=True)

print('| Fix | Item | Predicted | Debug measured | Interval | Rule |\n|---|---|---|---|---|---|'); print('\n'.join(rows))
print('\n### Debug vs Release beyond tolerance\n\n| Value | Release (review) | Debug | Δ | Tolerance |\n|---|---|---|---|---|'); print('\n'.join(flags) if flags else '| (none) | | | | |')
print('\n### Lag check (|lag| ≤ 48 smp required by the pre-registration, else reported)\n')
print(', '.join(f'{t} {q["lag"]:+d}' for t, q in sorted(M.items()) if 'lag' in q and abs(q['lag']) > 48))
