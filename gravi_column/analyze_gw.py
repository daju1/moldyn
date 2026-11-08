#!/usr/bin/env python3
"""Разбор логов gravi_walls (v1 и v2).

  python3 analyze_gw.py LOG [LOG ...] [--skip PS] [--calib CALIB_LOG] [--plot]

--skip PS   отбросить начало (установление), по умолчанию 200 пс
--calib     лог калибровочного прогона (стенки с перепадом), из него kappa
--plot      сохранить графики LOG.png (нужен matplotlib)

Поток тепла: наклон Q(t) после --skip, ошибка — по 5 блокам.
q > 0 : стенка отдаёт тепло газу.  q_up = (q_lo - q_hi)/2 : поток вверх через столб.
"""
import sys, numpy as np
NA = 6.0221367e23; R = 8.314510
A = (5e-9) * (5e-9)                      # площадь стенки, м^2 (ячейка 5x5 нм)
M = 0.0201797                            # кг/моль, Ne
CV = 1.5 * R; CP = 2.5 * R

def load(fn):
    hdr = [l for l in open(fn) if l.startswith('#')]
    a = np.array([l.split() for l in open(fn) if not l.startswith('#') and l.strip()], float)
    v2 = any('v2' in h for h in hdr)
    g = 0.0
    for h in hdr:
        for tok in h.split():
            if tok.startswith('g='): g = float(tok[2:])
    ncol = 11 if v2 else 7
    nb = (a.shape[1] - ncol) // 3
    prof = a[:, ncol:].reshape(len(a), nb, 3)
    return dict(a=a, v2=v2, g=g, nb=nb, prof=prof, t=a[:, 1])

def rate(t, Q, nblk=5):
    p = np.polyfit(t, Q, 1)[0]
    idx = np.array_split(np.arange(len(t)), nblk)
    r = [np.polyfit(t[i], Q[i], 1)[0] for i in idx if len(i) > 2]
    err = np.std(r, ddof=1) / np.sqrt(len(r)) if len(r) > 1 else np.nan
    return p, err

def analyze(fn, skip):
    d = load(fn); a, t = d['a'], d['t']
    m = t >= skip
    if m.sum() < 10: print(f'== {fn}: мало данных после {skip} пс (есть до {t[-1]:.1f} пс)'); return None
    conv = 1e3 / NA / A / 1e-12            # кДж/моль на систему за пс -> Вт/м^2
    qlo, elo = rate(t[m], a[m, 2]); qhi, ehi = rate(t[m], a[m, 3])
    qlo *= conv; qhi *= conv; elo *= conv; ehi *= conv
    prof = d['prof'][m]
    if d['v2']:
        T = prof[:, :, 2].mean(0); n = prof[:, :, 1].mean(0)
    else:                                   # v1: накопленные средние — берём разность последних
        T = prof[-1, :, 2]; n = prof[-1, :, 1]
    z = prof[0, :, 0]
    inner = slice(2, -2)
    k, kerr = np.polyfit(z[inner], T[inner], 1, cov=True); slope, s_err = k[0], np.sqrt(kerr[0, 0])
    gSI = abs(d['g']) * 1e12
    res = dict(fn=fn, q_lo=qlo, e_lo=elo, q_hi=qhi, e_hi=ehi, q_up=(qlo - qhi) / 2,
               e_up=0.5 * np.hypot(elo, ehi), slope=slope, s_err=s_err, g=gSI, T=T, n=n, z=z, d=d, m=m)
    print(f'== {fn}  (v2={d["v2"]}, g={gSI:.2e} м/с^2, окно {t[m][0]:.0f}..{t[m][-1]:.0f} пс)')
    print(f'  q_lo (нижняя отдаёт газу) = {qlo:+.3e} ± {elo:.1e} Вт/м^2')
    print(f'  q_hi (верхняя отдаёт газу)= {qhi:+.3e} ± {ehi:.1e} Вт/м^2')
    print(f'  поток вверх через столб    = {res["q_up"]:+.3e} ± {res["e_up"]:.1e} Вт/м^2')
    print(f'  dT/dz в объёме             = {slope:+.2f} ± {s_err:.2f} K/нм')
    if gSI > 0:
        print(f'     ориентиры: Циолковский -Mg/Cv = {-M*gSI/CV*1e-9:+.2f}, адиабата -Mg/Cp = {-M*gSI/CP*1e-9:+.2f}, изотерма 0 K/нм')
    print('  T(z): ' + ' '.join(f'{x:.0f}' for x in T))
    print('  n(z): ' + ' '.join(f'{x:.1f}' for x in n))
    if d['v2']:
        dr = a[m, 10]
        print(f'  дрейф энергии (Ekin+Epot+Egrav-E0-Q): {dr[0]:+.2f} .. {dr[-1]:+.2f} кДж/моль; Ekin ср = {a[m,7].mean():.1f}')
    return res

def main():
    args = sys.argv[1:]; skip = 200.0; calib = None; plot = False; logs = []
    i = 0
    while i < len(args):
        if args[i] == '--skip': skip = float(args[i+1]); i += 2
        elif args[i] == '--calib': calib = args[i+1]; i += 2
        elif args[i] == '--plot': plot = True; i += 1
        else: logs.append(args[i]); i += 1
    kappa = None
    if calib:
        c = analyze(calib, skip)
        if c is not None and abs(c['slope']) > 3 * c['s_err']:
            kappa = -c['q_up'] / (c['slope'] * 1e9)
            print(f'  -> kappa = -q_up/(dT/dz) = {kappa:.3e} Вт/(м K)')
    out = [analyze(f, skip) for f in logs]
    for r in out:
        if r is None or kappa is None or r['g'] == 0: continue
        gt = r['slope'] * 1e9
        q_four = -kappa * gt
        q_user = -kappa * (gt + M * r['g'] / CV)
        print(f'== сравнение для {r["fn"]}:')
        print(f'  измерено  q_up = {r["q_up"]:+.3e} ± {r["e_up"]:.1e} Вт/м^2')
        print(f'  Фурье     q = -k dT/dz            = {q_four:+.3e}')
        print(f'  ноутбук   q = -k (dT/dz + Mg/Cv)  = {q_user:+.3e}')
    if plot:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        for r in out:
            if r is None: continue
            a, t = r['d']['a'], r['d']['t']
            fig, ax = plt.subplots(1, 3, figsize=(14, 4))
            ax[0].plot(t, a[:, 2], label='Q_lo'); ax[0].plot(t, a[:, 3], label='Q_hi')
            ax[0].axvline(skip, ls=':', c='k'); ax[0].set_xlabel('t, пс'); ax[0].set_ylabel('кДж/моль'); ax[0].legend()
            ax[1].plot(r['z'], r['T'], 'o-'); ax[1].set_xlabel('z, нм'); ax[1].set_ylabel('T, K')
            if r['g'] > 0:
                zz = r['z']; T0 = r['T'][len(zz)//2]
                ax[1].plot(zz, T0 - M*r['g']/CV*1e-9*zz, '--', label='Циолковский')
                ax[1].plot(zz, T0 - M*r['g']/CP*1e-9*zz, ':', label='адиабата'); ax[1].legend()
            ax[2].plot(r['z'], r['n'], 'o-'); ax[2].set_xlabel('z, нм'); ax[2].set_ylabel('N в слое')
            fig.tight_layout(); fig.savefig(r['fn'] + '.png', dpi=110); print('  график:', r['fn'] + '.png')

if __name__ == '__main__': main()
