Design of [OpenMP]-based *dynamically scheduled Barrier-free* **Dynamic**
[PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)][page]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)][pagerank]. This is similar to barrierless non-blocking
implementations of the PageRank algorithm by Hemalatha Eedi et al. [(3)][eedi].
As ranks are updated in the same vector (with each iteration), the order in
which vertices are processed *affects* the final result (hence the adjective
*ordered*). However, as PageRank is an iteratively converging algorithm, results
obtained with either approach are *mostly the same*. **Barrier-free PageRank**
is an *ordered* *PageRank* where each thread processes a subset of vertices in
the graph independently, *without* waiting (with a barrier) for other threads to
complete an iteration. This minimizes unnecessary waits and allows each thread
to be on a *different iteration number* (which may or may not be beneficial for
convergence) [(3)][eedi].

Dynamic graphs, which change with time, have many applications. Computing ranks
of vertices from scratch on every update (*static PageRank*) may not be good
enough for an *interactive system*. In such cases, we only want to process ranks
of vertices which are likely to have changed. To handle any new vertices
added/removed, we first *adjust* the *previous ranks* (before the graph
update/batch) with a *scaled 1/N-fill* approach [(4)][adjust-ranks]. Then, with
**naive** **dynamic approach** we simply run the PageRank algorithm with the
*initial ranks* set to the adjusted ranks. Alternatively, with the (fully)
**dynamic approach** we first obtain a *subset of vertices* in the graph which
are likely to be affected by the update (using BFS/DFS from changed vertices),
and then perform PageRank computation on *only* this *subset of vertices*.

For each experiment below, we take fixed graphs as input, and generate random
batches of size `10^-8 |E|` to `0.1 |E|`, consisting of an equal mix of edge
deletions and insertions. For each batch size, we generate five different
(random) batches for averaging. A *schedule* of `dynamic, 2048` is used for
*OpenMP-based PageRank* as obtained in [(5)][pagerank-openmp]. We use the
follwing PageRank parameters: damping factor `α = 0.85`, tolerance `τ = 10^-10`,
and limit the maximum number of iterations to `L = 500.` The error between the
current and the previous iteration is obtained with *Li-norm*, and is used to
detect convergence. *Dead ends* in the graph are handled by adding self-loops to
all vertices in the graph (*loopall* approach [(6)][teleport]). Error in ranks
obtained for each approach is measured relative to the *sequential static
approach* using *Li-norm*.

The input data used for this experiment is available from the
[SuiteSparse Matrix Collection]. This experiment was done with guidance from
[Prof. Kishore Kothapalli], [Prof. Hemalatha Eedi], and [Prof. Sathya Peri].

[page]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[eedi]: https://ieeexplore.ieee.org/document/9407114
[teleport]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[adjust-ranks]: https://gist.github.com/wolfram77/eb7a3b2e44e3c2069e046389b45ead03
[pagerank]: https://github.com/puzzlef/pagerank
[pagerank-openmp]: https://github.com/puzzlef/pagerank-openmp-adjust-schedule
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Hemalatha Eedi]: https://jntuhceh.ac.in/faculty_details/5/dept/369
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis

<br>


### In the Absence of Faults

In this experiment, we compare the performance of *Static*, *Naive-dynamic*, and
*Dynamic Frontier* based *With-barrier* and *Barrier-free* PageRank. From the
results, we observe that our proposed *Dynamic Frontier* based *Barrier-free*
PageRank is `6.5x`, `5.8x`, `4.1x`, `2.4x`, and `1.4x` faster than the
*Naive-dynamic* approach for batch fractions of `10^-8` to `10^-4`. However,
beyond a batch size of `10^-3 |E|` (`1,000,000` edges for a billion-edge graph),
the performance drops below *Naive-dynamic*, and *Static* approaches (for both
*With-barrier* and *Barrier-free* PageRank). Such large batch sizes result in
*nearly* all of the vertices getting marked as *affected*, and hence the
performance drops. See [charts][charts1] included below, generated from
[sheets][sheets1].

[![](https://i.imgur.com/KpWMej4.png)][sheetp1]

[charts1]: https://imgur.com/a/hvq7z6q
[sheets1]: https://docs.google.com/spreadsheets/d/1Gjzv9drtd_zqOYYD_PqvSde_mcmQcnAf_CFg3FHbNwY/edit?usp=sharing
[sheetp1]: https://docs.google.com/spreadsheets/d/e/2PACX-1vSaT6GPh5wpb6ytZV9QKK1k3vjqq_x1rpuMTPd2k8b_GsNaVMKPXii4ATIHTU9BOqD8mxYV-i2RW143/pubhtml

<br>


### Strong Scaling

In this experiment, we study the strong-scaling behavior of *Dynamic Frontier*
based *With-barrier* and *Dynamic Frontier* based *Barrier-free* PageRank on
batch updates of a fixed size of `10^-4 |E|` in the absence of faults. Here, we
measure the speedup of each algorithm with an increasing number of threads from
`1` to `32` in multiples of `2` with respect to a single-threaded execution of
the algorithm. We additionally compare *Static* and *Naive-dynamic* versions of
*With-barrier* and *Barrier-free* PageRank.

From the results, we observe that all algorithms offer similar scaling
performance. When using `32` threads, *Dynamic Frontier* based *Barrier-free*
PageRank offers a average speedup of `15\times`, *Naive-dynamic* and *Static*
*Barrier-free* PageRank offer a average speedup of `14\times` and approaches
based on *With-barrier* PageRank offer an average speedup of `13.4\times`. This
demonstrates that *Dynamic Frontier* based *Barrier-free* PageRank offers good
scaling performance, that continues to be in line with that provided by the
PageRank algorithm in general. See [charts][charts6] included below, generated
from [sheets][sheets6].

[![](https://i.imgur.com/NkSxOJF.png)][sheetp6]
[![](https://i.imgur.com/JfAbMLR.png)][sheetp6]
[![](https://i.imgur.com/KVXOnpV.png)][sheetp6]

[charts6]: https://imgur.com/a/B8gwsnS
[sheets6]: https://docs.google.com/spreadsheets/d/1hgRbPQXj_O__m5nBcClo5IkJTM5vztklcCMSB6zU4yA/edit?usp=sharing
[sheetp6]: https://docs.google.com/spreadsheets/d/e/2PACX-1vS270_kwg8bQoapizI1vBRQzdD9eEiyDBaYKzz1tQzTo3XruQuRpDx_lT553pnt9GlScnGlC-ceSdi9/pubhtml

<br>


### With Uniform Thread delays

In this experiment, we compare the performance of *Dynamic Frontier* based
*With-barrier* and *Barrier-free* PageRank in the presence of random thread
delays that are uniformly distributed among all threads. From the results, we
observe that *Dynamic Frontier* based *With-barrier* PageRank is significantly
affected by an increasing delay probability. In contrast, *Dynamic Frontier*
based *Barrier-free* PageRank is relatively less affected. At a delay
probability of `10^-6`, *Dynamic Frontier* based *Barrier-free* PageRank
continues to be on average `1.5x`, `1.8x`, and `2.2x` faster on delay durations
of `50`, `100`, and `200` ms, respectively. See [charts][charts2] included
below, generated from [sheets][sheets2].

[![](https://i.imgur.com/MkvMiXK.png)][sheetp2]

[charts2]: https://imgur.com/a/GMu3xSw
[sheets2]: https://docs.google.com/spreadsheets/d/1YtzLia-sNlK9mJCtAnIHV5yG_sV1zSZUUloqQaHa9YQ/edit?usp=sharing
[sheetp2]: https://docs.google.com/spreadsheets/d/e/2PACX-1vT887D9hyBJRNYL9sW8ZgavC7pLjEoFde0KRjB4J7x_92wN3caMwgyGGrzvYjE6FtzBycBKVeSnziS3/pubhtml

<br>


### With Non-uniform Thread crashes

In this experiment, we measure the performance of *Dynamic Frontier* based
*Barrier-free PageRank* in the presence of random thread crashes that are
non-uniformly distributed (i.e., only the first `T` threads are affected by
random thread crashes) with a high probability of `10^-6`. From the results, we
observe that the performance of *Dynamic Frontier* based *Barrier-free* PageRank
drops at an exponential rate with an increase in the number of crashed threads.
When `28` out of `32` threads have crashed, it is, on average, `4.5x` slower and
fails to complete only when all `32` threads crash. See [charts][charts5]
included below, generated from [sheets][sheets5].

[![](https://i.imgur.com/vcr3OTE.png)][sheetp5]

[charts5]: https://imgur.com/a/vINN31x
[sheets5]: https://docs.google.com/spreadsheets/d/1h2qMsFpu81fDzIcLrNvDJkmjfjum1otZ6Oyn39CtvEI/edit?usp=sharing
[sheetp5]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTS1-_CFClCa9Ouvl080QnJ8dcLgyA_RlhdOMj-VRgwP7SwnZhk4nSBgpZArv2VctmbahTu_FfSTvM7/pubhtml

<br>


### Dynamic Contracting Frontier approach

We also experiment with a [Dynamic Contracting Frontier] approach, where we
remove vertices from the affected set once their ranks are computed. However,
this showed a small performance drop compared to the *Dynamic Frontier*
approach, and thus we do not discuss it any further.

[Dynamic Contracting Frontier]: https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/approach-cfrontier

<br>


### With Check and mark

With [Check and mark], we mark a vertex as affected only if its not already
marked, i.e., we check before marking. It is a very simple optimization, and we
apply it to *Dynamic Frontier* based *With-barrier* and *Barrier-free* PageRank.
Results indicate that this offer a small performance improvement (for large
batch updates).

[Check and mark]: https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/with-check-and-mark

<br>


### With Chunked mark

With [Chunked mark], we consider vertices to be grouped by IDs of size `2^n`,
and mark entire chunk as affected instead of marking just a single vertex. This
reduces the memory occupied by the affected flag vector, but has the side effect
of marking additional (unaffected) vertices are affected. It is a very simple
optimization, and we apply it to *Dynamic Frontier* based *With-barrier* and
*Barrier-free* PageRank. Results indicate that this does not offer *any*
improvement.

[Chunked mark]: https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/with-chunked-mark

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [Why Segmentation fault is happening in this openmp code?](https://stackoverflow.com/a/13266595/1413259)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)
- [Assigning default values to shell variables with a single command in bash](https://stackoverflow.com/a/28085062/1413259)
- [How to define a string literal in gcc command line?](https://stackoverflow.com/a/15220280/1413259)
- [How to split strings over multiple lines in Bash?](https://stackoverflow.com/a/46806081/1413259)
- [OMP_STACKSIZE :: OPENMP API Specification: Version 5.0 November 2018](https://www.openmp.org/spec-html/5.0/openmpse54.html)
- [Macros with values :: Linuxtopia](https://www.linuxtopia.org/online_books/an_introduction_to_gcc/gccintro_35.html)

<br>
<br>


[![](https://img.youtube.com/vi/OP-uxSvHUn8/maxresdefault.jpg)](https://www.youtube.com/watch?v=OP-uxSvHUn8)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/532937318.svg)](https://zenodo.org/badge/latestdoi/532937318)
