# Default settings: checkouts in $HOME, output in data/ at the repository root.
# Override any of these from the environment, or point CHEM_MACHINE at a file of
# your own that sets the same variables.
AK=${ATHENAK_ROOT:-$HOME/athenak}
AP=${ATHENAPP_ROOT:-$HOME/athena}
AKBIN=${AKBIN:-$AK/build/src/athena}
APBIN=${APBIN:-$AP/bin/athena}
DATA=${CHEMDATA:-$HERE/../data}
PY=${PYATHENA:-python3}
NP=${NP:-8}
