#include <utility>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
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




void runPagerankBatch(const string& data, size_t batch, size_t skip, int repeat) {
  using K = int;
  using T = TYPE;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
  vector<T> ranksOld, ranksAdj;
  vector<T> *initStatic  = nullptr;
  vector<T> *initDynamic = &ranksAdj;
  auto fl = [](auto u) { return true; };
  float damping   = 0.85;
  float tolerance = 1e-8;

  OutDiGraph<K> xo;
  stringstream  stream(data);
  while (true) {
    // Lets skip some edges.
    if (!readSnapTemporalW(xo, stream, skip)) break;
    auto x  = selfLoop(xo, None(), fl);
    auto xt = transposeWithDegree(x);
    auto a0 = pagerankMonolithicSeq<false, true>(x, xt, initStatic, {1, Li, damping, tolerance});
    auto ksOld = vertexKeys(x);
    ranksOld   = a0.ranks;

    // Read batch to be processed.
    auto yo = duplicate(xo);
    if (!readSnapTemporalW(yo, stream, batch)) break;
    auto y  = selfLoop(yo, None(), fl);
    auto yt = transposeWithDegree(y);
    auto a1 = pagerankMonolithicSeq<false, true>(y, yt, initStatic, {1, Li, damping, tolerance});
    auto ks = vertexKeys(y);

    do {
      // Find pagerank accelerated with OpenMP (static, unordered, no dead ends).
      auto a2 = pagerankMonolithicOmp<false, false>(y, yt, initStatic, {repeat, Li, damping, tolerance});
      auto e2 = l1Norm(a2.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnorderedStatic\n", y.order(), y.size(), a2.time, a2.iterations, e2);
      // Find pagerank accelerated with OpenMP (static, ordered, no dead ends).
      auto a3 = pagerankMonolithicOmp<true, false>(y, yt, initStatic, {repeat, Li, damping, tolerance});
      auto e3 = l1Norm(a3.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrderedStatic\n", y.order(), y.size(), a3.time, a3.iterations, e3);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (static, ordered, no dead ends, full error).
      auto a4 = pagerankBarrierfreeOmp<true, false, true>(y, yt, initStatic, {repeat, Li, damping, tolerance});
      auto e4 = l1Norm(a4.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmpStatic\n", y.order(), y.size(), a4.time, a4.iterations, e4);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (static, ordered, no dead ends, partial error).
      auto a5 = pagerankBarrierfreeOmp<true, false, false>(y, yt, initStatic, {repeat, Li, damping, tolerance});
      auto e5 = l1Norm(a5.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmpStatic\n", y.order(), y.size(), a5.time, a5.iterations, e5);
    } while (0);

    // Adjust ranks for dynamic Pagerank.
    ranksAdj.resize(y.span());
    adjustRanks(ranksAdj, ranksOld, ksOld, ks, 0.0f, float(ksOld.size())/ks.size(), 1.0f/ks.size());

    do {
      // Find pagerank accelerated with OpenMP (naive dynamic, unordered, no dead ends).
      auto a2 = pagerankMonolithicOmp<false, false>(y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e2 = l1Norm(a2.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnorderedNaiveDynamic\n", y.order(), y.size(), a2.time, a2.iterations, e2);
      // Find pagerank accelerated with OpenMP (naive dynamic, ordered, no dead ends).
      auto a3 = pagerankMonolithicOmp<true, false>(y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e3 = l1Norm(a3.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrderedNaiveDynamic\n", y.order(), y.size(), a3.time, a3.iterations, e3);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (naive dynamic, ordered, no dead ends, full error).
      auto a4 = pagerankBarrierfreeOmp<true, false, true>(y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e4 = l1Norm(a4.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmpNaiveDynamic\n", y.order(), y.size(), a4.time, a4.iterations, e4);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (naive dynamic, ordered, no dead ends, partial error).
      auto a5 = pagerankBarrierfreeOmp<true, false, false>(y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e5 = l1Norm(a5.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmpNaiveDynamic\n", y.order(), y.size(), a5.time, a5.iterations, e5);
    } while (0);

    do {
      // Find pagerank accelerated with OpenMP (dynamic, unordered, no dead ends).
      auto a2 = pagerankMonolithicOmpDynamic<false, false>(x, xt, y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e2 = l1Norm(a2.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpUnorderedDynamic\n", y.order(), y.size(), a2.time, a2.iterations, e2);
      // Find pagerank accelerated with OpenMP (dynamic, ordered, no dead ends).
      auto a3 = pagerankMonolithicOmpDynamic<true, false>(x, xt, y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e3 = l1Norm(a3.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrderedDynamic\n", y.order(), y.size(), a3.time, a3.iterations, e3);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (dynamic, ordered, no dead ends, full error).
      auto a4 = pagerankBarrierfreeOmpDynamic<true, false, true>(x, xt, y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e4 = l1Norm(a4.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeFullOmpDynamic\n", y.order(), y.size(), a4.time, a4.iterations, e4);
      // Find pagerank with barrier-free iterations accelerated with OpenMP (dynamic, ordered, no dead ends, partial error).
      auto a5 = pagerankBarrierfreeOmpDynamic<true, false, false>(x, xt, y, yt, initDynamic, {repeat, Li, damping, tolerance});
      auto e5 = l1Norm(a5.ranks, a1.ranks);
      printf("[%zu order; %zu size; %09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreePartOmpDynamic\n", y.order(), y.size(), a5.time, a5.iterations, e5);
    } while (0);

    // Now time to move on to next batch.
    xo = move(yo);
  }
}


void runPagerank(const string& data, int repeat) {
  size_t M = countLines(data), steps = 10;
  printf("Temporal edges: %zu\n\n", M);
  for (size_t batch=100; batch<=1000000; batch*=10) {
    size_t skip = max(int64_t(M/steps) - int64_t(batch), 0L);
    printf("# Batch size %.0e\n", double(batch));
    runPagerankBatch(data, batch, skip, repeat);
    printf("\n");
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  printf("Using graph %s ...\n", file);
  string data = readFile(file);
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(data, repeat);
  printf("\n");
  return 0;
}
