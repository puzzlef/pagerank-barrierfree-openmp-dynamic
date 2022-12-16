#!/usr/bin/env bash
src="pagerank-barrierfree-pthread"
out="$HOME/Logs/$src.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
rm -rf $src
git clone https://github.com/ionicf/$src
cd $src

# Run
g++ -std=c++17 -O3 -fopenmp -pthread main.cxx
stdbuf --output=L ./a.out ~/Data/web-Stanford.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-BerkStan.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-Google.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-NotreDame.mtx     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Slashdot0811.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Slashdot0902.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Epinions1.mtx     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coAuthorsDBLP.mtx     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coAuthorsCiteseer.mtx 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-LiveJournal1.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coPapersCiteseer.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coPapersDBLP.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/indochina-2004.mtx    2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/italy_osm.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/great-britain_osm.mtx 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/germany_osm.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/asia_osm.mtx          2>&1 | tee -a "$out"
