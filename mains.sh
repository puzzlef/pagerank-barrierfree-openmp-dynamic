#!/usr/bin/env bash
src="pagerank-barrierfree-openmp-static-vs-dynamic"
ulimit -s unlimited

# Download program
rm -rf $src
git clone https://github.com/puzzlef/$src
cd $src

# 1. Static vs Dynamic Barrier-free PageRank
export MAX_THREADS="24"
./main.sh

# 2. With uniform sleep failure
export FAILURE_TYPE="sleep"
export FAILURE_THREADS_BEGIN="$MAX_THREADS"
export FAILURE_THREADS_END="$MAX_THREADS"
./main.sh "--uniform-sleep"

# 3. With uniform crash failure
export FAILURE_TYPE="crash"
export FAILURE_THREADS_BEGIN="$MAX_THREADS"
export FAILURE_THREADS_END="$MAX_THREADS"
./main.sh "--uniform-crash"

# 4. With non-uniform sleep failure
export FAILURE_TYPE="sleep"
export FAILURE_THREADS_BEGIN="0"
export FAILURE_THREADS_END="$MAX_THREADS"
./main.sh "--nonuniform-sleep"

# 4. With non-uniform crash failure
export FAILURE_TYPE="crash"
export FAILURE_THREADS_BEGIN="0"
export FAILURE_THREADS_END="$MAX_THREADS"
./main.sh "--nonuniform-crash"
