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


### With Check and mark

With Check and mark, we mark a vertex as affected only if its not already
marked, i.e., we check before marking. It is a very simple optimization, and we
apply it to *Dynamic Frontier* based *With-barrier* and *Barrier-free* PageRank.
Results indicate that this offer a small performance improvement (for large
batch updates). Output is listed in [gist], and results are generated from
[sheets].

[![](https://i.imgur.com/6hkPDCh.png)][sheetp]
[![](https://i.imgur.com/BBF7Vt9.png)][sheetp]
[![](https://i.imgur.com/iIUgQqo.png)][sheetp]

[gist]: https://gist.github.com/wolfram77/a9a03a085ff959ae41868ff9017ce95b
[sheets]: https://docs.google.com/spreadsheets/d/1QtjYzbb5HTZ-I1_Sb18pAc7nAICILb4b0BKkvwc2HQs/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTpS5XjwdHEoptfZRwGL5juXDLCSLJQf55cjGo1ET6ZmqrBXlKyxyhQW8b-b1GqKUUULUO7XSiWf0w-/pubhtml

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
