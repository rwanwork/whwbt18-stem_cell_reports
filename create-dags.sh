#!/bin/bash

CORES=20
DAGS_DIR="Experiments/dags/"

TMP_FILE=`mktemp -p .`
grep -v combo7 Snakefile.py >${TMP_FILE}
snakemake --snakefile ${TMP_FILE} --cores ${CORES} --dag -np | dot -Tsvg >${DAGS_DIR}/combo7.svg
rm -f ${TMP_FILE}

