Design of [OpenMP]-based **dynamically scheduled Barrier-free** [PageRank algorithm]
for [link analysis].

For each experiment below, we take fixed graphs as input, and generate random
batches of size `10^-6 |E|` to `0.1 |E|`, consisting of edge insertions. For
each batch size, we generate five different (random) batches for averaging. A
*schedule* of `dynamic, 2048` is used for *OpenMP-based PageRank* as obtained in
[(1)][pagerank-openmp]. We use the follwing PageRank parameters: damping factor
`α = 0.85`, tolerance `τ = 10^-6`, and limit the maximum number of iterations to
`L = 500.` The error between the current and the previous iteration is obtained
with *L1-norm*, and is used to detect convergence. *Dead ends* in the graph are
handled by adding self-loops to all vertices in the graph (*loopall* approach
[(2)][teleport]). Error in ranks obtained for each approach is measured relative
to the *sequential static approach* using *L1-norm*.

The input data used for this experiment is available from the
[SuiteSparse Matrix Collection]. This experiment was done with guidance from
[Prof. Kishore Kothapalli], [Prof. Hemalatha Eedi], and [Prof. Sathya Peri].

[teleport]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[pagerank-openmp]: https://github.com/puzzlef/pagerank-openmp
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
*Dynamic* (*Traversal*, *DFS-based*, DFS time *not* included) *Barrier-free*
*PageRank*. From the results, we observe that **dynamically scheduled**
**Barrier-free PageRank** is **faster** than *With-Barrier* PageRank. See [charts][charts1]
included below, generated from [sheets][sheets1].

[![](https://i.imgur.com/s5oct3o.png)][sheetp1]

[charts1]: https://imgur.com/a/SnknAwf
[sheets1]: https://docs.google.com/spreadsheets/d/1HrefgcERVI-g3Q3aTQIUtP8T5zywaNdHlmeJB2jgsOM/edit?usp=sharing
[sheetp1]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTMuvH38QJ0w9rv_v1bZ1IJywEYCWcrekxa17UzoPsTjHOOAXJ20m8vidhX_Eua6StynTCsLcwaIwuk/pubhtml

<br>


### With Uniform Thread delays

In this experiment, we compare the *Static*, *Naive-dynamic*, and *Dynamic*
*Barrier-free PageRank* in the presence of random thread delays that are
uniformly distributed among all threads. From the results, we observe that
**dynamically scheduled Barrier-free PageRank** is significantly **less**
**impacted** by random thread delays than *With-Barrier* PageRank. See
[charts][charts2] included below, generated from [sheets][sheets2].

[![](https://i.imgur.com/9lmx6Xs.png)][sheetp2]
[![](https://i.imgur.com/phOhUFY.png)][sheetp2]
[![](https://i.imgur.com/fkiZQIF.png)][sheetp2]

[charts2]: https://imgur.com/a/BJ9zwju
[sheets2]: https://docs.google.com/spreadsheets/d/16XrJ7dgTwhL0O4XTrMDNCD0Ho_q5G_LWGxDc4G3yfFU/edit?usp=sharing
[sheetp2]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTkMCH2q2ONfBtmqxd4LdEzVwMs1GjuAX59DNPubzvprl8ETErWqauh03QLtlZB78E6PJk7mBYlDz4D/pubhtml

<br>


### With Non-uniform Thread delays

In this experiment, we compare the *Static*, *Naive-dynamic*, and *Dynamic*
*Barrier-free PageRank* in the presence of random thread delays that are
non-uniformly distributed (i.e., only the first `T` threads are affected by
random thread delays). From the results, we observe that **dynamically**
**scheduled Barrier-free PageRank** is significantly **less impacted** by
random thread delays than *With-Barrier* PageRank. See [charts][charts3]
included below, generated from [sheets][sheets3].

[![](https://i.imgur.com/FCUUm34.png)][sheetp3]
[![](https://i.imgur.com/k0jmpm5.png)][sheetp3]
[![](https://i.imgur.com/J5ocZ10.png)][sheetp3]
[![](https://i.imgur.com/U8OFOlB.png)][sheetp3]
[![](https://i.imgur.com/EB9pX9v.png)][sheetp3]
[![](https://i.imgur.com/MD99pEU.png)][sheetp3]
[![](https://i.imgur.com/uO0KevV.png)][sheetp3]

[charts3]: https://imgur.com/a/yP7OBvr
[sheets3]: https://docs.google.com/spreadsheets/d/1XueRFUavO0IA1PCFIUAYIfloWnJwPuMpAmz3xB2aJ7s/edit?usp=sharing
[sheetp3]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTr6kHOP2CWYBMFm0yBArkFHpOE3UZl4lYyqn6AfuNGkYxwLx5bkkcSQjpu2vFV4MHUSmdgiRk1OVO6/pubhtml

<br>


### With Uniform Thread crashes

In this experiment, we compare the *Static*, *Naive-dynamic*, and *Dynamic*
*Barrier-free PageRank* in the presence of random thread crashes that are
uniformly distributed among all threads. From the results, we observe that
**dynamically scheduled Barrier-free PageRank** is significantly **less**
**impacted** by random thread crashes, while *With-Barrier* PageRank fails
to complete PageRank computation (even when a single thread crashes). See
[charts][charts4] included below, generated from [sheets][sheets4].

[![](https://i.imgur.com/A1WTvGT.png)][sheetp4]

[charts4]: https://imgur.com/a/Ez19YW9
[sheets4]: https://docs.google.com/spreadsheets/d/1jdxgSaSz4SFr7yVp0F_LWrGJBtx6A0OEt3MlEIDsR7I/edit?usp=sharing
[sheetp4]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRsGwM2vzcXI3HjaNEjmMghyd2rvlYXNx7LuIoEFDvjOybsulgToiG_3Mx5HfKeekuBIyD42PtCA2tD/pubhtml

<br>


### With Non-uniform Thread crashes

In this experiment, we compare the *Static*, *Naive-dynamic*, and *Dynamic*
*Barrier-free PageRank* in the presence of random thread crashes that are
non-uniformly distributed (i.e., only the first `T` threads are affected by
random thread crashes). From the results, we observe that **dynamically**
**scheduled Barrier-free PageRank** is significantly **less impacted** by random
thread crashes, while *With-Barrier* PageRank fails to complete PageRank
computation (even when a single thread crashes). See [charts][charts5] included
below, generated from [sheets][sheets5].

[![](https://i.imgur.com/WjORP8z.png)][sheetp5]
[![](https://i.imgur.com/qdXLcpy.png)][sheetp5]
[![](https://i.imgur.com/3Uk5vui.png)][sheetp5]
[![](https://i.imgur.com/LcGRgeJ.png)][sheetp5]
[![](https://i.imgur.com/FeCXTM5.png)][sheetp5]
[![](https://i.imgur.com/AfsOHrl.png)][sheetp5]
[![](https://i.imgur.com/c2RAk4Q.png)][sheetp5]

[charts5]: https://imgur.com/a/GpzFK0Y
[sheets5]: https://docs.google.com/spreadsheets/d/1krcSYLZbW7OtIQrSA44bRpR4-UpuSkb6SQwmgWCdU9M/edit?usp=sharing
[sheetp5]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTqNPDs3RRzwsPPJW4uMyVnNWlkbXk9py1BdwQ2A6wSlLm6EzFVtxIj6jCFdByg9jL5WHT2LKcMPqgR/pubhtml

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
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/532937318.svg)](https://zenodo.org/badge/latestdoi/532937318)
