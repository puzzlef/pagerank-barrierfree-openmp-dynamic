#!/usr/bin/env bash
src="pagerank-barrierfree-openmp-static-vs-dynamic"
ulimit -s unlimited

# Download program
if [[ "$DOWNLOAD" != "0" ]]; then
  rm -rf $src
  git clone https://github.com/puzzlef/$src
  cd $src
fi

# Don't need to download program again.
export DOWNLOAD="0"

# 1. Static vs Dynamic Barrier-free PageRank
export MAX_THREADS="24"
./main.sh

# For uniform failure
export BATCH_INSERTIONS_BEGIN="0.001"
export BATCH_INSERTIONS_END="0.001"
export FAILURE_THREADS_BEGIN="$MAX_THREADS"
export FAILURE_THREADS_END="$MAX_THREADS"

# 2. With uniform sleep failure
export FAILURE_TYPE="sleep"
./main.sh "--uniform-sleep"

# 3. With uniform crash failure
export FAILURE_TYPE="crash"
./main.sh "--uniform-crash"

# For non-uniform failure
export BATCH_INSERTIONS_BEGIN="0.001"
export BATCH_INSERTIONS_END="0.001"
export FAILURE_DURATION_BEGIN="100"
export FAILURE_DURATION_END="100"
export FAILURE_THREADS_BEGIN="0"
export FAILURE_THREADS_END="$MAX_THREADS"

# 4. With non-uniform sleep failure
export FAILURE_TYPE="sleep"
./main.sh "--nonuniform-sleep"

# 5. With non-uniform crash failure
export FAILURE_TYPE="crash"
./main.sh "--nonuniform-crash"

# Signal completion
event="puzzlef_${src//-/_}"
curl -X POST https://maker.ifttt.com/trigger/${event}/json/with/key/${IFTTT_KEY}
