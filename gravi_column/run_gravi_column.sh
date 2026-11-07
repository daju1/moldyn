#!/bin/bash
# Столб неона в поле тяжести между двумя тепловыми стенками (patch gravi_walls).
# Запускать ВНУТРИ контейнера (podman_run.sh), из этой папки.
set -e -x

RELEASE_DIR=$(dirname $(dirname ${PWD}))/winlibghemical3/linux/release20111012

export LD_LIBRARY_PATH="${RELEASE_DIR}/local/lib"
HERE="$(cd "$(dirname "$0")" && pwd)"
EXE="${RELEASE_DIR}/examples/example_moldyn"

run_case () {   # $1=имя  $2=T_низ  $3=T_верх  $4=g
  export GW_ENABLE=1 GW_G=$4 GW_ZLO=-8.0 GW_ZHI=8.0 GW_TLO=$2 GW_THI=$3 GW_NBINS=16 GW_EVERY=5000
  export GW_OUT="${HERE}/$1.gw.log"
  # -E : без термостата в объёме (только стенки);  -B 0 : без оптимизации ячейки
  # -T $2 : начальное распределение Максвелла;  -t 2 fs;  5e6 шагов = 10 нс
  "${EXE}" -E -B 0 -T $2 -t 2.0 -N 5000000 -f 50000 -v 15 \
      "${HERE}/ne_column.gpr" "${HERE}/ne_column.box" "${HERE}/$1.traj" > "${HERE}/$1.stdout" 2>&1
}

case "${1:-equal}" in
  equal)  run_case equal_T300_g6   300 300 -6.0 ;;   # главный опыт: равные стенки, тяжесть
  nog)    run_case equal_T300_g0   300 300  0.0 ;;   # контроль: без тяжести (должен быть ноль)
  calib)  run_case calib_330_270_g0 330 270 0.0 ;;   # калибровка kappa: перепад 60 K без тяжести
  all)    "$0" equal & "$0" nog & "$0" calib & wait ;;
esac
