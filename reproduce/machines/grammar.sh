# grammar (KIAS CPU cluster): AthenaK CPU+MPI and Athena++ with CVODE.
# Shares $HOME and /gpfs with syntax. Run inside a Slurm allocation; bare mpirun
# on the login node is for the smoke test only.
#
# Not yet run here. Build first:
#   AthenaK:  ~/ai-notes/scripts/prepare_athenak.py --device cpu --mpi
#   Athena++: configure.py --prob=read_vtk --flux=hllc --chemistry=gow17
#             --cvode_path=<sundials> -mpi, with grammar_modules.sh sourced
#
# Submit from the repository root:
#   sbatch -p normal -N 1 -n 8 --mem-per-cpu=2G -t 01:00:00 \
#     -o /gpfs/jeonggyukim/athenak-chem/logs/turb64-%j.out reproduce/run_turb64.sh
AK=${ATHENAK_ROOT:-$HOME/athenak-sweep}
AP=${ATHENAPP_ROOT:-$HOME/athena-pp}
AKBIN=${AKBIN:-$AK/build-cpu/src/athena}
APBIN=${APBIN:-$AP/bin/athena_mpi}
DATA=${CHEMDATA:-/gpfs/jeonggyukim/athenak-chem/runs/turb64}
PY=${PYATHENA:-$HOME/miniconda3/envs/pyathena/bin/python}
NP=${NP:-8}
source "${AI_SCRIPTS:-$HOME/ai-notes/scripts}/grammar_modules.sh"
