#!/usr/bin/env bash
src="pagerank-barrierfree-openmp-static-vs-dynamic"
out="$HOME/Logs/$src.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
rm -rf $src
git clone https://github.com/puzzlef/$src
cd $src

# Fixed config
: "${TYPE:=double}"
: "${MAX_THREADS:=24}"
: "${REPEAT_BATCH:=5}"
: "${REPEAT_METHOD:=1}"
# Parameter sweep for batch (randomly generated)
: "${BATCH_UNIT:=\%}"
: "${BATCH_DELETIONS_BEGIN:=0}"
: "${BATCH_DELETIONS_END:=0}"
: "${BATCH_DELETIONS_STEP:=+=1}"
: "${BATCH_INSERTIONS_BEGIN:=0.000001}"
: "${BATCH_INSERTIONS_END:=0.001001}"
: "${BATCH_INSERTIONS_STEP:=\*=10}"
# Parameter sweep for failure (simulated)
: "${FAILURE_TYPE:=none}"
: "${FAILURE_DURATION_BEGIN:=1}"
: "${FAILURE_DURATION_END:=100}"
: "${FAILURE_DURATION_STEP:=\*=10}"
: "${FAILURE_PROBABILITY_BEGIN:=0.0001}"
: "${FAILURE_PROBABILITY_END:=1}"
: "${FAILURE_PROBABILITY_STEP:=\*=10}"
: "${FAILURE_THREADS_BEGIN:=0}"
: "${FAILURE_THREADS_END:=0}"
: "${FAILURE_THREADS_STEP:=+=1}"
# Define macros (dont forget to add here)
DEFINES=""\
"-DTYPE=$TYPE"\
"-DMAX_THREADS=$MAX_THREADS"\
"-DREPEAT_BATCH=$REPEAT_BATCH"\
"-DREPEAT_METHOD=$REPEAT_METHOD"\
"-DBATCH_TYPE=\"$BATCH_TYPE\""\
"-DBATCH_UNIT=\"$BATCH_UNIT\""\
"-DBATCH_DELETIONS_BEGIN=$BATCH_DELETIONS_BEGIN"\
"-DBATCH_DELETIONS_END=$BATCH_DELETIONS_END"\
"-DBATCH_DELETIONS_STEP=$BATCH_DELETIONS_STEP"\
"-DBATCH_INSERTIONS_BEGIN=$BATCH_INSERTIONS_BEGIN"\
"-DBATCH_INSERTIONS_END=$BATCH_INSERTIONS_END"\
"-DBATCH_INSERTIONS_STEP=$BATCH_INSERTIONS_STEP"\
"-DFAILURE_TYPE=\"$FAILURE_TYPE\""\
"-DFAILURE_DURATION_BEGIN=$FAILURE_DURATION_BEGIN"\
"-DFAILURE_DURATION_END=$FAILURE_DURATION_END"\
"-DFAILURE_DURATION_STEP=$FAILURE_DURATION_STEP"\
"-DFAILURE_PROBABILITY_BEGIN=$FAILURE_PROBABILITY_BEGIN"\
"-DFAILURE_PROBABILITY_END=$FAILURE_PROBABILITY_END"\
"-DFAILURE_PROBABILITY_STEP=$FAILURE_PROBABILITY_STEP"\
"-DFAILURE_THREADS_BEGIN=$FAILURE_THREADS_BEGIN"\
"-DFAILURE_THREADS_END=$FAILURE_THREADS_END"\
"-DFAILURE_THREADS_STEP=$FAILURE_THREADS_STEP"

# Run
g++ "$DEFINES" -std=c++17 -O3 -fopenmp main.cxx
stdbuf --output=L ./a.out ~/Data/indochina-2004.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2002.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/arabic-2005.mtx     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2005.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/webbase-2001.mtx    2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/it-2004.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/sk-2005.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-LiveJournal.mtx 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Orkut.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/asia_osm.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/europe_osm.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_A2a.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_V1r.mtx        2>&1 | tee -a "$out"
