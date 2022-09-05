#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class H>
void runPagerank(const G& x, const H& xt, int repeat) {
  using T = TYPE;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
  vector<T> *init = nullptr;
  float damping   = 0.85;

  // Use L1-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find pagerank using a single thread for reference (unordered, no dead ends).
    auto a0 = pagerankMonolithicSeq<false, false>(x, xt, init, {1, L1, damping, tolerance});
    // Find pagerank accelerated with OpenMP (unordered, no dead ends).
    auto a1 = pagerankMonolithicOmp<false, false>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered       {tol_norm: L1, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find pagerank accelerated with OpenMP (ordered, no dead ends).
    auto a2 = pagerankMonolithicOmp<true, false>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered         {tol_norm: L1, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, full error).
    auto a3 = pagerankBarrierfreeOmp<true, false, true>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, partial error).
    auto a4 = pagerankBarrierfreeOmp<true, false, false>(x, xt, init, {repeat, L1, damping, tolerance});
    auto e4 = l1Norm(a4.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: %.0e}\n", a4.time, a4.iterations, e4, tolerance);
  }

  // Use L2-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find pagerank using a single thread for reference (unordered, no dead ends).
    auto a0 = pagerankMonolithicSeq<false, false>(x, xt, init, {1, L2, damping, tolerance});
    // Find pagerank accelerated with OpenMP (unordered, no dead ends).
    auto a1 = pagerankMonolithicOmp<false, false>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered       {tol_norm: L2, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find pagerank accelerated with OpenMP (ordered, no dead ends).
    auto a2 = pagerankMonolithicOmp<true, false>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered         {tol_norm: L2, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, full error).
    auto a3 = pagerankBarrierfreeOmp<true, false, true>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmp {tol_norm: L2, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, partial error).
    auto a4 = pagerankBarrierfreeOmp<true, false, false>(x, xt, init, {repeat, L2, damping, tolerance});
    auto e4 = l1Norm(a4.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmp {tol_norm: L2, tolerance: %.0e}\n", a4.time, a4.iterations, e4, tolerance);
  }

  // Use Li-norm for convergence check.
  for (float tolerance=1e-1; tolerance>=1e-15; tolerance/=10) {
    // Find pagerank using a single thread for reference (unordered, no dead ends).
    auto a0 = pagerankMonolithicSeq<false, false>(x, xt, init, {1, Li, damping, tolerance});
    // Find pagerank accelerated with OpenMP (unordered, no dead ends).
    auto a1 = pagerankMonolithicOmp<false, false>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e1 = l1Norm(a1.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnordered       {tol_norm: Li, tolerance: %.0e}\n", a1.time, a1.iterations, e1, tolerance);
    // Find pagerank accelerated with OpenMP (ordered, no dead ends).
    auto a2 = pagerankMonolithicOmp<true, false>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e2 = l1Norm(a2.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered         {tol_norm: Li, tolerance: %.0e}\n", a2.time, a2.iterations, e2, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, full error).
    auto a3 = pagerankBarrierfreeOmp<true, false, true>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e3 = l1Norm(a3.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmp {tol_norm: Li, tolerance: %.0e}\n", a3.time, a3.iterations, e3, tolerance);
    // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends, partial error).
    auto a4 = pagerankBarrierfreeOmp<true, false, false>(x, xt, init, {repeat, Li, damping, tolerance});
    auto e4 = l1Norm(a4.ranks, a0.ranks);
    printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmp {tol_norm: Li, tolerance: %.0e}\n", a4.time, a4.iterations, e4, tolerance);
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  printf("Loading graph %s ...\n", file);
  auto x  = readMtxOutDiGraph(file); println(x);
  auto fl = [](auto u) { return true; };
  selfLoopU(x, None(), fl); print(x); printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegree(x);  print(xt); printf(" (transposeWithDegree)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(x, xt, repeat);
  printf("\n");
  return 0;
}
