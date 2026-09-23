#!/bin/bash
# Decaying 64^3 turbulence with GOW17 chemistry: AthenaK semi-implicit against
# Athena++ CVODE, both started from one shared VTK snapshot.
#
#   bash reproduce/run_turb64.sh                 every stage, in order
#   bash reproduce/run_turb64.sh si cvode join   only the named stages
#   SMOKE=1 bash reproduce/run_turb64.sh         every stage, 20 cycles each
#
# Stages:
#   drive  AthenaK driven turbulence, hydro only, t = 0 to 12 (1 rank, ~5 min)
#   ic     last drive64 dump -> shared_ic64.vtk, Athena 4.2 legacy VTK
#   si     AthenaK, semi-implicit chemistry, t = 0 to 2.0454 (NP ranks, ~40 s)
#   bdf    AthenaK, Kokkos BDF at the CVODE tolerances, same interval (NP ranks)
#   cvode  Athena++, CVODE chemistry, same interval (NP ranks, ~18 min)
#   join   merge the Athena++ per-MeshBlock VTK files of each dump into joined/
#   plot   history from the .hst files, history from the dumps, and slice
#          comparisons (final frame and movie)
# Times are for 8 ranks on the MacBook.
#
# Machine settings come from machines/<machine>.sh, picked by hostname or by
# MACHINE=mac|grammar|syntax. Every setting there can be overridden from the
# environment (ATHENAK_ROOT, ATHENAPP_ROOT, AKBIN, APBIN, CHEMDATA, PYATHENA, NP).
#
# On the KIAS clusters submit from the repository root, since Slurm runs a
# spooled copy of this file and cannot locate reproduce/ from it:
#   grammar: sbatch -p normal -N 1 -n 8 --mem-per-cpu=2G -t 01:00:00 \
#              -o /gpfs/jeonggyukim/athenak-chem/logs/turb64-%j.out \
#              reproduce/run_turb64.sh
#   syntax:  sbatch -p h200 --gres=gpu:1 -n 1 -t 00:30:00 \
#              -o /gpfs/jeonggyukim/athenak-chem/logs/turb64-%j.out \
#              reproduce/run_turb64.sh drive ic si
#
# Output goes to $CHEMDATA. A stage refuses to overwrite a run directory that
# already holds run.out; move it aside or point CHEMDATA elsewhere.
set -e

if [ -n "$SLURM_JOB_ID" ]; then
  HERE=$SLURM_SUBMIT_DIR/reproduce
else
  HERE=$(cd "$(dirname "$0")" && pwd)
fi

if [ -z "$MACHINE" ]; then
  case $(hostname -s) in
    grammar*) MACHINE=grammar ;;
    syntax|syn[0-9]*) MACHINE=syntax ;;
    *MacBook*) MACHINE=mac ;;
    *) echo "ERROR  unknown host $(hostname -s); set MACHINE=mac|grammar|syntax" >&2
       exit 1 ;;
  esac
fi
source "$HERE/machines/$MACHINE.sh"

# The smoke test runs every stage for 20 cycles, into its own directory, with
# figures beside the data rather than in the repository.
OVR=()
FIGDIR=${FIGDIR:-$HERE/../figures}
if [ -n "$SMOKE" ]; then
  [ -n "$CHEMDATA" ] || DATA=$DATA-smoke
  OVR=(time/nlim=20)
  FIGDIR=$DATA/figures
fi

STAGES=${*:-drive ic si bdf cvode join plot}
want() { case " $STAGES " in *" $1 "*) return 0 ;; esac; return 1; }

mkdir -p "$DATA"

# A figure is reproducible only if the code that made its numbers is recorded.
# Both working trees carry uncommitted changes (the read_vtk pgens), so the
# commit alone is not enough: save the diff and the untracked sources too.
stamp() {
  local name=$1 repo=$2
  {
    echo "$name: $(git -C "$repo" rev-parse --short HEAD)" \
         "on $(git -C "$repo" rev-parse --abbrev-ref HEAD)"
    git -C "$repo" status --porcelain | sed 's/^/  /'
  } >> "$DATA/versions.txt"
  git -C "$repo" diff HEAD > "$DATA/$name.diff"
  git -C "$repo" ls-files --others --exclude-standard -- src > "$DATA/$name-untracked.txt"
  if [ -s "$DATA/$name-untracked.txt" ]; then
    tar -C "$repo" -cf "$DATA/$name-untracked.tar" -T "$DATA/$name-untracked.txt"
  fi
}

check_fresh() {
  local newer
  newer=$(find "$AK/src" -newer "$AKBIN" \( -name '*.cpp' -o -name '*.hpp' \) | head -1)
  if [ -n "$newer" ]; then
    echo "ERROR  $newer is newer than $AKBIN; rebuild first." >&2
    exit 1
  fi
}

check_ap() {
  local cfg
  cfg=$("$APBIN" -c)
  if ! echo "$cfg" | grep -qE 'Riemann solver: +hllc'; then
    echo "ERROR  $APBIN is not built with --flux=hllc; AthenaK uses hllc." >&2
    exit 1
  fi
  if ! echo "$cfg" | grep -qE 'Problem generator: +read_vtk'; then
    echo "ERROR  $APBIN is not built with --prob=read_vtk." >&2
    exit 1
  fi
  echo "$cfg" > "$DATA/athena-pp-configure.txt"
}

prepare() {
  local dir=$DATA/$1
  if [ -e "$dir/run.out" ]; then
    echo "ERROR  $dir already holds a run; move it aside first." >&2
    exit 1
  fi
  mkdir -p "$dir"
  cp -f "$HERE/inputs/$2" "$dir/"
}

run() {
  local dir=$DATA/$1 np=$2
  shift 2
  echo "=== $1 in $dir ($np ranks)"
  ( cd "$dir" && { time mpirun -np "$np" "$@" > run.out 2> run.err ; } 2> time.txt )
  tail -3 "$dir/time.txt"
}

echo "date: $(date '+%Y-%m-%d %H:%M'), machine: $MACHINE, stages: $STAGES" \
  "${SMOKE:+(smoke)}" >> "$DATA/versions.txt"
stamp athenak "$AK"
stamp athena-pp "$AP"

if want drive; then
  check_fresh
  prepare drive64_run drive64.athinput
  run drive64_run 1 "$AKBIN" -i drive64.athinput "${OVR[@]}"
fi

if want ic; then
  last=$(ls "$DATA"/drive64_run/bin/drive64.hydro_w.*.bin | tail -1)
  "$PY" "$HERE/athenak_bin_to_vtk.py" "$last" "$DATA/shared_ic64.vtk"
  echo "shared_ic64.vtk from $last" >> "$DATA/versions.txt"
fi

if want si; then
  check_fresh
  prepare turb64_si_run turb64_si.athinput
  cp -f "$DATA/shared_ic64.vtk" "$DATA/turb64_si_run/"
  run turb64_si_run "$NP" "$AKBIN" -i turb64_si.athinput "${OVR[@]}"
fi

if want bdf; then
  check_fresh
  prepare turb64_bdf_run turb64_bdf.athinput
  cp -f "$DATA/shared_ic64.vtk" "$DATA/turb64_bdf_run/"
  run turb64_bdf_run "$NP" "$AKBIN" -i turb64_bdf.athinput "${OVR[@]}"
fi

if want cvode; then
  check_ap
  prepare turb64_cvode_run turb64_cvode.athinput
  cp -f "$DATA/shared_ic64.vtk" "$DATA/turb64_cvode_run/"
  run turb64_cvode_run "$NP" "$APBIN" -i turb64_cvode.athinput "${OVR[@]}"
fi

if want join; then
  dir=$DATA/turb64_cvode_run
  nblock=$(ls "$dir"/turb64_cvode.block*.out1.00000.vtk | wc -l)
  mkdir -p "$dir/joined"
  for f in "$dir"/turb64_cvode.block0.out1.*.vtk; do
    n=${f##*.out1.}
    blocks=()
    for (( b = 0; b < nblock; b++ )); do blocks+=("$dir/turb64_cvode.block$b.out1.$n"); done
    "$AP/vis/vtk/join_vtk++" -o "$dir/joined/turb64_cvode.block0.out1.$n" "${blocks[@]}" > /dev/null
  done
  echo "joined $(ls "$dir/joined" | wc -l | tr -d ' ') dumps of $((nblock)) MeshBlocks each"
fi

if want plot; then
  mkdir -p "$FIGDIR"
  "$PY" "$HERE/../scripts/plot_turb_hst.py" \
    --ak "$DATA/turb64_si_run/turb64_si.hydro.hst" \
    --ap "$DATA/turb64_cvode_run/turb64_cvode.hst" \
    --out "$FIGDIR/fig_turb_hst.png"
  "$PY" "$HERE/../scripts/plot_turb_history_dumps.py" \
    --ak-si-dir "$DATA/turb64_si_run" --ak-si-base turb64_si \
    --ak-bdf-dir "$DATA/turb64_bdf_run" --ak-bdf-base turb64_bdf \
    --ap-cvode-dir "$DATA/turb64_cvode_run" --ap-cvode-base turb64_cvode \
    --out "$FIGDIR/fig_turb_history_dumps.png"
  "$PY" "$HERE/../scripts/plot_turb_evolution.py" \
    --ak-dir "$DATA/turb64_si_run" --ak-base turb64_si \
    --ap-dir "$DATA/turb64_cvode_run/joined" --ap-base turb64_cvode \
    --outdir "$FIGDIR"
fi
