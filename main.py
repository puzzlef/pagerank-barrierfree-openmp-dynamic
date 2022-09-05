# https://www.kaggle.com/wolfram77/puzzlef-pagerank-barrierfree-openmp-static-vs-dynamic
import os
from IPython.display import FileLink
src="pagerank-barrierfree-openmp-static-vs-dynamic"
inp="/kaggle/input/graphs"
out="{}.txt".format(src)
!printf "" > "$out"
display(FileLink(out))
!echo ""

# Download program
!rm -rf $src
!git clone https://github.com/puzzlef/$src

# Run
!g++ -std=c++17 -O3 -fopenmp $src/main.cxx
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/email-Eu-core-temporal.txt 2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/CollegeMsg.txt             2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/sx-mathoverflow.txt        2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/sx-askubuntu.txt           2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/sx-superuser.txt           2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/wiki-talk-temporal.txt     2>&1 | tee -a "$out"
!ulimit -s unlimited && stdbuf --output=L ./a.out ~/data/sx-stackoverflow.txt       2>&1 | tee -a "$out"
