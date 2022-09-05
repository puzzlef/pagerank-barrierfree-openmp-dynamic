Effect of using different values of tolerance with barrier-free iterations in
[OpenMP]-based ordered [PageRank algorithm] for [link analysis].

`TODO`

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)]. This is similar to barrierless non-blocking implementations
of the PageRank algorithm by Hemalatha Eedi et al. [(3)]. As ranks are updated
in the same vector (with each iteration), the order in which vertices are
processed *affects* the final result (hence the adjective *ordered*). However,
as PageRank is an iteratively converging algorithm, results obtained with either
approach are *mostly the same*. **Barrier-free PageRank** is an *ordered*
*PageRank* where each thread processes a subset of vertices in the graph
independently, *without* waiting (with a barrier) for other threads to complete an
iteration. This minimizes unnecessary waits and allows each thread to be on a
*different iteration number* (which may or may not be beneficial for convergence)
[(3)].

In this experiment, we perform two diffrent approaches of barrier-free
iterations of *OpenMP-based ordered PageRank*; one in which each thread detects
convergence by measuring the difference between the previous and the current
ranks of all the vertices (**full**), and the other in which the difference is
measured between the previous and current ranks of only the subset of vertices
being processed by each thread (**part**). This is done while adjusting the
tolerance `τ` from `10^-1` to `10^-14` with three different tolerance functions:
`L1-norm`, `L2-norm`, and `L∞-norm`. We also compare it with OpenMP-based
unordered and ordered PageRank for the same tolerance and tolerance function. We
use a damping factor of `α = 0.85` and limit the maximum number of iterations to
`L = 500`. The error between the approaches is calculated with *L1-norm*. The
*sequential unordered* approach is considered to be the *gold standard* (wrt to
which error is measured). *Dead ends* in the graph are handled by always
teleporting any vertex in the graph at random (*teleport* approach [(4)]). The
teleport contribution to all vertices is calculated *once* (for all vertices) at
the begining of each iteration.

From the results, we observe that the **partial-error based barrier-free**
**PageRank** is **faster than full-error version**, but generally has **similar**
**performance as that of standard OpenMP-based PageRank (unordered, ordered)**,
although it **requires the fewest iterations to converge**. However the,
**partial barrier-free approach** also **generally yields greater error in**
**ranks** (with respect to sequential unordered approach).

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli],
[Prof. Dip Sankar Banerjee], and [Prof. Sathya Peri].

<br>

```bash
$ g++ -std=c++17 -O3 -fopenmp main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 2594400 [directed] {} (selfLoopAllVertices)
# order: 281903 size: 2594400 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00004.980 ms; 004 iters.] [0.0000e+00 err.] pagerankOmpUnordered       {tol_norm: L1, tolerance: 1e-01}
# [00007.029 ms; 004 iters.] [7.2516e-02 err.] pagerankOmpOrdered         {tol_norm: L1, tolerance: 1e-01}
# [00007.962 ms; 003 iters.] [6.1728e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-01}
# [00003.015 ms; 001 iters.] [3.6976e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-01}
# [00009.053 ms; 009 iters.] [0.0000e+00 err.] pagerankOmpUnordered       {tol_norm: L1, tolerance: 1e-02}
# [00012.224 ms; 008 iters.] [2.0893e-02 err.] pagerankOmpOrdered         {tol_norm: L1, tolerance: 1e-02}
# [00015.114 ms; 008 iters.] [2.1312e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-02}
# [00006.250 ms; 003 iters.] [1.0583e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-02}
# ...
# [00385.129 ms; 500 iters.] [0.0000e+00 err.] pagerankOmpUnordered       {tol_norm: Li, tolerance: 1e-14}
# [00132.396 ms; 098 iters.] [1.9295e-07 err.] pagerankOmpOrdered         {tol_norm: Li, tolerance: 1e-14}
# [00170.530 ms; 108 iters.] [1.9585e-07 err.] pagerankBarrierfreeFullOmp {tol_norm: Li, tolerance: 1e-14}
# [00123.138 ms; 094 iters.] [1.9663e-07 err.] pagerankBarrierfreePartOmp {tol_norm: Li, tolerance: 1e-14}
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 8285825 [directed] {} (selfLoopAllVertices)
# order: 685230 size: 8285825 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00008.499 ms; 004 iters.] [0.0000e+00 err.] pagerankOmpUnordered       {tol_norm: L1, tolerance: 1e-01}
# [00009.376 ms; 004 iters.] [1.0193e-01 err.] pagerankOmpOrdered         {tol_norm: L1, tolerance: 1e-01}
# [00013.724 ms; 004 iters.] [1.0509e-01 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-01}
# [00004.638 ms; 001 iters.] [3.4068e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-01}
# [00016.418 ms; 009 iters.] [0.0000e+00 err.] pagerankOmpUnordered       {tol_norm: L1, tolerance: 1e-02}
# [00017.874 ms; 008 iters.] [2.9763e-02 err.] pagerankOmpOrdered         {tol_norm: L1, tolerance: 1e-02}
# [00026.116 ms; 009 iters.] [2.8237e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-02}
# [00008.815 ms; 003 iters.] [1.0316e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-02}
# ...
```

[![](https://i.imgur.com/f2EYibQ.png)][sheetp]
[![](https://i.imgur.com/L2Pqr8R.png)][sheetp]
[![](https://i.imgur.com/lOJyiPG.png)][sheetp]
[![](https://i.imgur.com/vdMuchv.png)][sheetp]

[![](https://i.imgur.com/hJUfY4I.png)][sheetp]
[![](https://i.imgur.com/LT5Yt1v.png)][sheetp]
[![](https://i.imgur.com/Pouaonx.png)][sheetp]
[![](https://i.imgur.com/f0xSHPU.png)][sheetp]

[![](https://i.imgur.com/oXne2DD.png)][sheetp]
[![](https://i.imgur.com/4W2G1IV.png)][sheetp]
[![](https://i.imgur.com/JgpwWjL.png)][sheetp]
[![](https://i.imgur.com/lQ9u5u9.png)][sheetp]

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)

<br>
<br>


[![](https://i.imgur.com/6UYGStl.jpg)](https://www.youtube.com/watch?v=2k2orAjC4aw)<br>
[![DOI](https://zenodo.org/badge/532185741.svg)](https://zenodo.org/badge/latestdoi/532185741)


[(1)]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[(2)]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[(3)]: https://ieeexplore.ieee.org/document/9407114
[(4)]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/4e88175074a8b134cddcd76d5b1b3542
[charts]: https://imgur.com/a/PjMk8zX
[sheets]: https://docs.google.com/spreadsheets/d/1bIdZGAOHDXN3bCjasOPAjOZCpWaTQpgYlZYv6SKeRT4/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vSo6nUPtx3ahAFElFW2bLjFiDWHCQ3NHLv_PuFaLjR-u7FySxlq54tdNdPveBEWw06WHODtl8HqJzUd/pubhtml
