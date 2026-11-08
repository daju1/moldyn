#!/bin/bash
# Столб неона в поле тяжести между двумя тепловыми стенками (патч gravi_walls v2).
# Запускать ВНУТРИ контейнера/на Linux-машине, из этой папки.
#   ./run_gravi_column.sh equal 1      главный опыт, зерно 1
#   ./run_gravi_column.sh equal 2      он же, независимая траектория (второй компьютер)
#   ./run_gravi_column.sh nog 1        контроль без тяжести
#   ./run_gravi_column.sh calib 1      калибровка kappa: 330/270 K без тяжести
#   ./run_gravi_column.sh calibg 1     то же с тяжестью (kappa в поле)
set -e -x
RELEASE_DIR="${RELEASE_DIR:-$(cd "$(dirname "$0")/../../winlibghemical/linux/release20111012" && pwd)}"
export LD_LIBRARY_PATH="${RELEASE_DIR}/local/lib"
HERE="$(cd "$(dirname "$0")" && pwd)"
EXE="${RELEASE_DIR}/examples/example_moldyn"

DT=${DT:-0.5}            # фс; 2 фс неустойчиво
NSTEPS=${NSTEPS:-4000000} # 4e6 x 0.5 фс = 2 нс
EVERY=${EVERY:-5000}      # строка лога каждые 2.5 пс

run_case () {   # $1=имя $2=T_низ $3=T_верх $4=g $5=seed
  export GW_ENABLE=1 GW_G=$4 GW_ZLO=-8.0 GW_ZHI=8.0 GW_TLO=$2 GW_THI=$3 GW_NBINS=16 GW_EVERY=${EVERY} GW_SEED=$5
  export GW_OUT="${HERE}/$1_s$5.gw.log"
  T0=$(( ($2 + $3) / 2 ))
  # -E: без термостата вообще (с патчем v2 — и без стадий нагрева); -B 0: без оптимизации ячейки
  "${EXE}" -E -B 0 -T ${T0} -t ${DT} -N ${NSTEPS} -f 200000 -v 15 \
      "${HERE}/${GPR:-ne_column.gpr}" "${HERE}/ne_column.box" "${HERE}/$1_s$5.traj" > "${HERE}/$1_s$5.stdout" 2>&1
}

SEED=${2:-1}
case "${1:-equal}" in
  equal)  run_case equal_T300_g6 300 300 -6.0 ${SEED} ;;
  nog)    run_case nog_T300_g0   300 300  0.0 ${SEED} ;;
  calib)  run_case calib_330_270_g0 330 270 0.0 ${SEED} ;;
  calibg) run_case calib_330_270_g6 330 270 -6.0 ${SEED} ;;
  # плотный вариант: 1600 атомов (n≈4 нм^-3, λ≈0.8 нм, H/λ≈20)
  dense)       GPR=ne_column_dense.gpr run_case dense_T300_g6        300 300 -6.0 ${SEED} ;;
  dense_nog)   GPR=ne_column_dense.gpr run_case dense_nog_T300_g0    300 300  0.0 ${SEED} ;;
  dense_calib) GPR=ne_column_dense.gpr run_case dense_calib_330_270_g0 330 270 0.0 ${SEED} ;;
  *) echo "unknown case $1"; exit 1 ;;
esac
