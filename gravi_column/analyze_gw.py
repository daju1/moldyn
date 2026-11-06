#!/usr/bin/env python3
"""Разбор логов gravi_walls: поток тепла через стенки и профиль T(z).
Запуск:  python3 analyze_gw.py equal_T300_g6.gw.log [ещё логи ...]
Единицы: Q в кДж/моль (на моль атомов-«событий»), площадь стенки из box.
"""
import sys, numpy as np
NA=6.0221367e23
def load(fn):
    rows=[l.split() for l in open(fn) if not l.startswith('#')]
    a=np.array(rows,float); return a
for fn in sys.argv[1:]:
    a=load(fn); t=a[:,1]*1e-12               # с
    Qlo=a[:,2]*1e3/NA; Qhi=a[:,3]*1e3/NA      # Дж (на систему)
    half=len(t)//2                             # первую половину отбрасываем (установление)
    def rate(Q): p=np.polyfit(t[half:],Q[half:],1); return p[0]
    A=(2*2.5e-9)*(2*2.5e-9)                   # площадь 5x5 нм
    qlo=rate(Qlo)/A; qhi=rate(Qhi)/A          # Вт/м^2, + = стенка отдаёт газу
    nb=(a.shape[1]-7)//3; prof=a[-1,7:].reshape(nb,3)
    z,nz,T=prof[:,0],prof[:,1],prof[:,2]
    k=np.polyfit(z[2:-2],T[2:-2],1)[0]
    print(f'== {fn}')
    print(f'  тепло от нижней стенки  q_lo = {qlo:+.3e} Вт/м2   (+ = в газ)')
    print(f'  тепло от верхней стенки q_hi = {qhi:+.3e} Вт/м2')
    print(f'  поток вверх через столб (q_lo - q_hi)/2 = {(qlo-qhi)/2:+.3e} Вт/м2')
    print(f'  наклон T(z) в объёме = {k:+.3f} K/нм ; T: '+' '.join(f'{x:.0f}' for x in T))
