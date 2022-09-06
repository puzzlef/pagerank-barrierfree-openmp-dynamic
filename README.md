Performance of static vs dynamic barrier-free iterations in [OpenMP]-based
ordered [PageRank algorithm] for [link analysis].

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

Dynamic graphs, which change with time, have many applications. Computing ranks
of vertices from scratch on every update (*static PageRank*) may not be good
enough for an *interactive system*. In such cases, we only want to process ranks
of vertices which are likely to have changed. To handle any new vertices
added/removed, we first *adjust* the *previous ranks* (before the graph
update/batch) with a *scaled 1/N-fill* approach [(4)]. Then, with **naive**
**dynamic approach** we simply run the PageRank algorithm with the *initial ranks*
set to the adjusted ranks. Alternatively, with the (fully) **dynamic approach**
we first obtain a *subset of vertices* in the graph which are likely to be
affected by the update (using BFS/DFS from changed vertices), and then perform
PageRank computation on *only* this *subset of vertices*.

In this experiment, we compare the performance of **static**, **naive dynamic**,
and (fully) **barrier-free iterations in dynamic OpenMP-based ordered PageRank**
(along with similar *unordered/ordered OpenMP-based approaches*). We take
*temporal graphs* as input, and add edges to our in-memory graph in batches of
size `10^2 to 10^6`. However we do *not* perform this on every point on the
temporal graph, but *only* on *5 time samples* of the graph (5 samples are good
enough to obtain an average). At each time sample we load `B` edges (where *B*
is the batch size), and perform *static*, *naive dynamic*, and *dynamic*
PageRank. At each time sample, each approach is performed *5* *times* to obtain
an average time for that sample.  A *schedule* of `dynamic, 2048` is used for
*OpenMP-based PageRank* as obtained in [(5)]. We use the follwing PageRank
parameters: damping factor `α = 0.85`, tolerance `τ = 10^-6`, and limit the
maximum number of iterations to `L = 500.` The error between the current and the
previous iteration is obtained with *L1-norm*, and is used to detect
convergence. *Dead ends* in the graph are handled by adding self-loops to all
vertices in the graph (*loopall* approach [(6)]). Error in ranks obtained for
each approach is measured relative to the *sequential static approach* using
*L1-norm*.

From the results, we observe that **OpenMP-based unordered PageRank is faster**
**than ordered approach**, and **Barrier-free PageRank with partial error**
**measurement is faster than the one with full error measurement**. **Partial**
**barrier-free dynamic PageRank is the fastest approach for batch sizes below or**
**equal to** `10^3`, while **OpenMP-based unordered dynamic PageRank is the**
**fastest for larger batch sizes**. We also note that OpenMP-based dynamic and
naive dynamic have almost identical performance, while barrier-free dynamic
approaches are faster what naive dynamic approaches. In terms of the number of
iterations, we observe that parital barrier-free dynamic/naive dynamic
approaches requires the fewest iterations to converge for all batch sizes. We
also note that all dynamic and naive dynamic approaches require similar number
of iterations to converge.

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli],
[Prof. Dip Sankar Banerjee], and [Prof. Sathya Peri].

<br>

```bash
$ g++ -std=c++17 -O3 -fopenmp main.cxx
$ ./a.out ~/data/email-Eu-core-temporal.txt
$ ./a.out ~/data/CollegeMsg.txt
$ ...

# Using graph /home/subhajit/data/email-Eu-core-temporal.txt ...
# OMP_NUM_THREADS=12
# Temporal edges: 332335
#
# # Batch size 1e+02
# [751 order; 7703 size; 00000.619 ms; 063 iters.] [0.0000e+00 err.] pagerankOmpUnorderedStatic
# [751 order; 7703 size; 00000.859 ms; 075 iters.] [5.6016e-07 err.] pagerankOmpOrderedStatic
# [751 order; 7703 size; 00020.857 ms; 089 iters.] [6.1278e-07 err.] pagerankBarrierfreeFullOmpStatic
# [751 order; 7703 size; 00006.136 ms; 067 iters.] [2.1334e-06 err.] pagerankBarrierfreePartOmpStatic
# [751 order; 7703 size; 00000.620 ms; 057 iters.] [9.0126e-07 err.] pagerankOmpUnorderedNaiveDynamic
# [751 order; 7703 size; 00000.675 ms; 057 iters.] [8.4747e-07 err.] pagerankOmpOrderedNaiveDynamic
# [751 order; 7703 size; 00006.907 ms; 062 iters.] [8.4526e-07 err.] pagerankBarrierfreeFullOmpNaiveDynamic
# [751 order; 7703 size; 00003.605 ms; 037 iters.] [2.2513e-06 err.] pagerankBarrierfreePartOmpNaiveDynamic
# [751 order; 7703 size; 00000.612 ms; 057 iters.] [9.0126e-07 err.] pagerankOmpUnorderedDynamic
# [751 order; 7703 size; 00000.659 ms; 057 iters.] [8.4745e-07 err.] pagerankOmpOrderedDynamic
# [751 order; 7703 size; 00005.945 ms; 062 iters.] [8.3778e-07 err.] pagerankBarrierfreeFullOmpDynamic
# [751 order; 7703 size; 00003.136 ms; 037 iters.] [2.2675e-06 err.] pagerankBarrierfreePartOmpDynamic
# ...
# [986 order; 25915 size; 00001.915 ms; 064 iters.] [0.0000e+00 err.] pagerankOmpUnorderedStatic
# [986 order; 25915 size; 00002.375 ms; 072 iters.] [5.7090e-07 err.] pagerankOmpOrderedStatic
# [986 order; 25915 size; 00027.283 ms; 080 iters.] [1.0225e-06 err.] pagerankBarrierfreeFullOmpStatic
# [986 order; 25915 size; 00020.180 ms; 064 iters.] [2.0462e-06 err.] pagerankBarrierfreePartOmpStatic
# [986 order; 25915 size; 00000.861 ms; 026 iters.] [1.3071e-06 err.] pagerankOmpUnorderedNaiveDynamic
# [986 order; 25915 size; 00000.957 ms; 026 iters.] [1.2276e-06 err.] pagerankOmpOrderedNaiveDynamic
# [986 order; 25915 size; 00010.578 ms; 027 iters.] [1.2453e-06 err.] pagerankBarrierfreeFullOmpNaiveDynamic
# [986 order; 25915 size; 00004.964 ms; 017 iters.] [1.8024e-06 err.] pagerankBarrierfreePartOmpNaiveDynamic
# [986 order; 25915 size; 00000.856 ms; 026 iters.] [1.3071e-06 err.] pagerankOmpUnorderedDynamic
# [986 order; 25915 size; 00000.955 ms; 026 iters.] [1.2276e-06 err.] pagerankOmpOrderedDynamic
# [986 order; 25915 size; 00007.328 ms; 028 iters.] [1.2456e-06 err.] pagerankBarrierfreeFullOmpDynamic
# [986 order; 25915 size; 00003.954 ms; 017 iters.] [1.8240e-06 err.] pagerankBarrierfreePartOmpDynamic
#
# # Batch size 1e+03
# [751 order; 7703 size; 00000.623 ms; 063 iters.] [0.0000e+00 err.] pagerankOmpUnorderedStatic
# [751 order; 7703 size; 00000.815 ms; 075 iters.] [5.6016e-07 err.] pagerankOmpOrderedStatic
# [751 order; 7703 size; 00012.586 ms; 091 iters.] [6.4177e-07 err.] pagerankBarrierfreeFullOmpStatic
# [751 order; 7703 size; 00008.811 ms; 067 iters.] [2.1239e-06 err.] pagerankBarrierfreePartOmpStatic
# [751 order; 7703 size; 00000.654 ms; 060 iters.] [1.0349e-06 err.] pagerankOmpUnorderedNaiveDynamic
# [751 order; 7703 size; 00000.707 ms; 060 iters.] [9.6145e-07 err.] pagerankOmpOrderedNaiveDynamic
# [751 order; 7703 size; 00011.240 ms; 073 iters.] [8.0149e-07 err.] pagerankBarrierfreeFullOmpNaiveDynamic
# [751 order; 7703 size; 00006.967 ms; 050 iters.] [1.9816e-06 err.] pagerankBarrierfreePartOmpNaiveDynamic
# [751 order; 7703 size; 00000.662 ms; 060 iters.] [1.0360e-06 err.] pagerankOmpUnorderedDynamic
# [751 order; 7703 size; 00000.715 ms; 060 iters.] [9.6246e-07 err.] pagerankOmpOrderedDynamic
# [751 order; 7703 size; 00010.556 ms; 068 iters.] [8.4651e-07 err.] pagerankBarrierfreeFullOmpDynamic
# [751 order; 7703 size; 00006.814 ms; 050 iters.] [2.0133e-06 err.] pagerankBarrierfreePartOmpDynamic
# ...
```

[![](https://i.imgur.com/wrcIxBc.png)][sheetp]
[![](https://i.imgur.com/Uz6lMb8.png)][sheetp]

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


[![](https://i.imgur.com/7Cuj7c9.jpg)](https://www.youtube.com/watch?v=OP-uxSvHUn8)<br>
[![DOI](https://zenodo.org/badge/532937318.svg)](https://zenodo.org/badge/latestdoi/532937318)


[(1)]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[(2)]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[(3)]: https://ieeexplore.ieee.org/document/9407114
[(4)]: https://gist.github.com/wolfram77/eb7a3b2e44e3c2069e046389b45ead03
[(5)]: https://github.com/puzzlef/pagerank-openmp-adjust-schedule
[(6)]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/7630b0484e14e9776d7d9fc68ed1db75
[charts]: https://imgur.com/a/7lnBmO9
[sheets]: https://docs.google.com/spreadsheets/d/1ETJfHz6i_zX7kdNTs6wCrr_qMzeduejr0ZQY8HPZ04c/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vQzhScoS8sxyx8uFXzIVO7qtt2cn2Vsb_yZg9ihpFcJQFxWS_Z1HBlgCR1snF18iuDYok3K5Pndx1Ct/pubhtml
