# Jeong-Gyus-MacBook-Pro: AthenaK CPU+MPI and Athena++ with CVODE.
AK=${ATHENAK_ROOT:-$HOME/Projects/athenak-chem}
AP=${ATHENAPP_ROOT:-$HOME/Projects/athena-pp}
AKBIN=${AKBIN:-$AK/build-cpu-mpi/src/athena}
APBIN=${APBIN:-$AP/bin/athena_mpi}
DATA=${CHEMDATA:-$HOME/Documents/Athenak-chem-tests/turb64}
PY=${PYATHENA:-/opt/homebrew/Caskroom/miniforge/base/envs/pyathena/bin/python}
NP=${NP:-8}
