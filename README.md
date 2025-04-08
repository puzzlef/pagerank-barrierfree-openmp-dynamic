Design of [OpenMP]-based *Barrier-free* **Dynamic** [PageRank algorithm] for
[link analysis].

**Barrier-free PageRank** [(1)][eedi] is an *asynchronous* PageRank variant
[(2)][pagerank] where each thread independently processes a subset of vertices
in the graph *without waiting* for other threads to complete an iteration. This
*minimizes unnecessary waits* and allows threads to be on different iteration
numbers. We add self-loops to each vertex in the graph to avoid the cost of
computing common teleport contribution to all vertices [(3)][teleport].

When dealing with *dynamic graphs* that change over time, computing ranks from
scratch on every update (static PageRank) may not be suitable for interactive
systems. In such cases, it is preferable to process only the ranks of vertices
likely to have changed. Here, we propose the **Dynamic Frontier approach** for
**With-barrier** and **Barrier-free PageRank** that *incrementally identifies an*
*active subset of vertices* in the graph which are likely to be affected by the
update, and then performs PageRank computation on only this subset of vertices.
This approach performs better than *Naive-dynamic* and *Traversal-based*
*(BFS/DFS) approaches*. We observe Traversal-based approaches to be slower in
general compared to Naive-dynamic for all but the smallest of batch sizes.
Accordingly, we do not include Traversal-based approach in our experiments.

The input data used for our experiments is available from the
[SuiteSparse Matrix Collection]. This experiment was done with guidance from
[Prof. Kishore Kothapalli], [Prof. Hemalatha Eedi], and [Prof. Sathya Peri].

[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Hemalatha Eedi]: https://jntuhceh.ac.in/faculty_details/5/dept/369
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[page]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[eedi]: https://ieeexplore.ieee.org/document/9407114
[teleport]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[adjust-ranks]: https://gist.github.com/wolfram77/eb7a3b2e44e3c2069e046389b45ead03
[pagerank]: https://github.com/puzzlef/pagerank
[pagerank-openmp]: https://github.com/puzzlef/pagerank-openmp-adjust-schedule

<br>


### In the Absence of Faults

In this experiment we compare **With-barrier** and **Barrier-free PageRank**
(*Static*, *Naive-dynamic*, and **Dynamic Frontier**) with a batch size of `B` =
`10^-8` to `0.1 |E|`. The batch update consists of a random batch of edge
deletions (of size `B`), and then we insert the same edges back. We compute
Dynamic PageRanks (*Naive-dynamic*, *Dynamic Frontier*) back to back (i.e., once
upon deletions, and then upon insertions), and then compare obtained ranks
wrt **reference ranks** on the original graph. Reference ranks are obtained
using a **very low tolerance** (`10^-100`, limited to `500` iterations)
*With-barrier* PageRank.

From results we observe that **Dynamic Frontier outperforms** Naive-dynamic upto
a batch size of `10^-3 |E|`. For a billion-edge graph, this amounts to edge
deletions of 1 million vertices, followed by edge insertions of 1 million
vertices. We also observe that **Barrier-free outperforms** With-barrier for all
batch sizes, which can be attributed to the lack of waiting time spent waiting
at the barrier. Further, we observe that **Barrier-free** obtains ranks
of **equivalent quality** as *With-barrier* (slightly better in case of
*Naive-dynamic*, *Dynamic Frontier*).

> See
> [code](https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/input-large),
> [output](https://gist.github.com/wolfram77/d78d33aa13b77156eea7c4c2bc2960e9), or
> [sheets][sheets-e1].

[![](https://i.imgur.com/VHol26s.png)][sheets-e1]
[![](https://i.imgur.com/2vByekl.png)][sheets-e1]

[sheets-e1]: https://docs.google.com/spreadsheets/d/1Gjzv9drtd_zqOYYD_PqvSde_mcmQcnAf_CFg3FHbNwY/edit?usp=sharing

<br>


### With Uniform Thread delays

In this experiment, we compare **With-barrier** and **Barrier-free PageRank**
(*Static*, *Naive-dynamic*, and **Dynamic Frontier**) in the presence of *random*
*thread sleep* on a batch size of `10^-4 |E|`. During rank computation of each
vertex, **a thread may randomly sleep** with a **probability** of `10^-9` to
`10^-6`. On a graph with 10 million vertices, this amounts to an average of
`0.01` to `10` threads sleeps per iteration. We try three different
sleep **durations**, i.e., `50ms`, `100ms`, `200ms` (JVM's GC pause is set to a
soft limit of [200ms][jvm-gc-pause] as reference). Batch updates are generated
as above.

Results indicate that **Barrier-free outperforms** *With-barrier* by a
significant margin when thread sleeps are more common. This can be useful in
non-dedicated systems where we occasionally run other activities in addition to
the main computation. Further, we observe that **Barrier-free** obtains ranks
of **equivalent quality** (slightly better) as *With-barrier*.

> See
> [code](https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/input-large),
> [output](https://gist.github.com/wolfram77/46d76338b3fca8f7661fcc7044eb8f20), or
> [sheets][sheets-e2].

[![](https://i.imgur.com/HyTeNWf.png)][sheets-e2]
[![](https://i.imgur.com/SkWUpdE.png)][sheets-e2]

[sheets-e2]: https://docs.google.com/spreadsheets/d/1YtzLia-sNlK9mJCtAnIHV5yG_sV1zSZUUloqQaHa9YQ/edit?usp=sharing
[jvm-gc-pause]: https://dzone.com/articles/how-to-reduce-long-gc-pause

<br>


### With Non-uniform Thread crashes

In this experiment, we test **Dynamic Frontier based Barrier-free PageRank** in
the presence of *random thread crash* on a batch size of `10^-4 |E|`. During
rank computation of each vertex, **a thread may randomly crash** (*crash-stop*
*model*) with a **high probability** of `10^-5`. On a graph with 10 million
vertices, this amounts to an average of `100` threads crashing per iteration. We
*limit* the number of *threads that may crash* from `0` to `56` threads. Batch
updates are generated as above. Note that *With-barrier PageRank cannot tolerate*
*thread crashes*.

Results indicate that **Barrier-free tolerates thread crashes** to the maximum
extent possible, with **almost no drop is result quality**. Note that runtime
taken by the algorithm increases with the number of crashed threads, which is
expected as the workload is distributed among the remaining threads. Further,
there is **no added cost** to offering this safety to crash-stop faults. We
anticipate this to be useful in robotic systems which operate in harsh
enviroments.

> See
> [code](https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/input-large),
> [output](https://gist.github.com/wolfram77/ca5fc1fb63ddeea4664186c0e89b4dc5), or
> [sheets][sheets-e3].

[![](https://i.imgur.com/ViqaUju.png)][sheets-e3]
[![](https://i.imgur.com/mHuwyFI.png)][sheets-e3]

[sheets-e3]: https://docs.google.com/spreadsheets/d/1f2w_uPXbKYyv5Mau_P0zgmv2zvWUXEI7gsT4CjJTa6w/edit?usp=sharing

<br>


### Strong scaling (in the absence of faults)

In this experiment, we test the strong-scaling behavior of **Dynamic Frontier**
**based Barrier-free PageRank** in the absence of faults on a batch size of
`10^-4 |E|`. We adjust the number of threads from `1` to `64` threads in multiples of `2`. Batch
updates are generated as above.

Results indicate that **Dynamic Frontier based Barrier-free PageRank** offers a
speedup of `21.3x` on `64` threads over `1` thread. We observe that the rate of
increase in speedup drops after `32` threads, which is likely due to the AMD CPU
being internally partitioned into 4 NUMA nodes.

> See
> [code](https://github.com/puzzlef/pagerank-barrierfree-openmp-dynamic/tree/input-large),
> [output](https://gist.github.com/wolfram77/48c714d9fbbd8a85372d2b2e1590dc19), or
> [sheets][sheets-e3].

[![](https://i.imgur.com/CIncReC.png)][sheets-e3]
[![](https://i.imgur.com/Npk3Ykp.png)][sheets-e3]

[sheets-e3]: https://docs.google.com/spreadsheets/d/1jhMEvnoHBifuZfRo5Qcy11A-N3xni7h0tVhFJBHNE30/edit?usp=sharing

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
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/pagerank-barrierfree-openmp-dynamic)
