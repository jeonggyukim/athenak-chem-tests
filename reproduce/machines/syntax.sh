# syntax (KIAS GPU cluster): AthenaK on one GPU. Athena++ CVODE belongs on
# grammar; syntax runs the drive, ic and si stages, and the output lands on the
# shared /gpfs where grammar can read it.
#
# Not yet run here. Build first, matching KOKKOS_ARCH to the partition
# (HOPPER90 for h100/h200, AMPERE80 for a100, AMPERE86 for a10/a40):
#   ~/ai-notes/scripts/prepare_athenak.py --device gpu --arch HOPPER90
#
# Submit from the repository root:
#   sbatch -p h200 --gres=gpu:1 -n 1 -t 00:30:00 \
#     -o /gpfs/jeonggyukim/athenak-chem/logs/turb64-%j.out \
#     reproduce/run_turb64.sh drive ic si
AK=${ATHENAK_ROOT:-$HOME/athenak-sweep}
AP=${ATHENAPP_ROOT:-$HOME/athena-pp}
AKBIN=${AKBIN:-$AK/build-gpu-${KOKKOS_ARCH:-HOPPER90}/src/athena}
APBIN=${APBIN:-$AP/bin/athena_mpi}
DATA=${CHEMDATA:-/gpfs/jeonggyukim/athenak-chem/runs/turb64}
PY=${PYATHENA:-$HOME/miniconda3/envs/pyathena/bin/python}
NP=${NP:-1}
set +u
if ! type module >/dev/null 2>&1; then
  export MODULEPATH=${MODULEPATH:-/opt/ohpc/pub/modulefiles}
  . /opt/ohpc/admin/lmod/lmod/init/bash >/dev/null
fi
module purge >/dev/null 2>&1
module load gnu12/12.2.0 openmpi4/4.1.5 cuda/12.6.3
