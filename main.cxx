#include <algorithm>
#include <chrono>
#include <random>
#include <thread>
#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE double
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 24
#endif




template <class G, class R>
auto addRandomEdges(G& a, R& rnd, size_t i, size_t n, size_t batchSize) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v);
    a.addEdge(v, u);
    insertions.push_back(make_tuple(u, v));
    insertions.push_back(make_tuple(v, u));
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(a, rnd, i, n, None(), fe); }, retries);
  updateOmpU(a);
  return insertions;
}


template <class G, class R>
auto removeRandomEdges(G& a, R& rnd, size_t i, size_t n, size_t batchSize) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    a.removeEdge(v, u);
    deletions.push_back(make_tuple(u, v));
    deletions.push_back(make_tuple(v, u));
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(a, rnd, i, n, fe); }, retries);
  updateOmpU(a);
  return deletions;
}




template <class G, class H>
void runExperiment(const G& x, const H& xt, int repeat) {
  using  K = typename G::key_type;
  using  V = TYPE;
  vector<V> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  size_t batchSize = size_t(0.001 * x.size());
  int    retries   = 5;
  // Get ranks of vertices on original graph (static).
  auto fu = [&](ThreadInfo *thread, auto v) {};
  auto a0 = pagerankBasicOmp(xt, init, {1}, fu);
  // Batch of additions only (dynamic).
  for (int sleepDuration=1; sleepDuration<=100; sleepDuration*=10) {
    for (float sleepProbability=0.0f; sleepProbability<0.101f; sleepProbability+=0.02f) {
      for (int batchCount=1; batchCount<=1; ++batchCount) {
        auto y  = duplicate(x);
        auto insertions = addRandomEdges(y, rnd, 1, x.span()-1, batchSize); vector<tuple<K, K>> deletions;
        auto yt = transposeWithDegreeOmp(y);
        auto fc = [](bool v) { return v==true; };
        size_t affectedCount = countIf(pagerankAffectedVerticesTraversal(x, deletions, insertions), fc);
        // Do something (sleep) after processing each vertex.
        chrono::milliseconds sd(sleepDuration);
        float sp = sleepProbability / x.order();
        auto  fv = [&](ThreadInfo *thread, auto v) {
          uniform_real_distribution<float> dis(0.0f, 1.0f);
          if (dis(thread->rnd) < sp) this_thread::sleep_for(sd);
        };
        // Follow a specific result logging format, which can be easily parsed later.
        auto flog  = [&](auto ans, auto ref, const char *technique) {
          auto ear = pagerankBasicOmp(yt, &ans.ranks, {1}, fu);
          auto err = l1NormOmp(ans.ranks, ref.ranks);
          LOG(
            "{+%.2e batch, %.2e aff, %03dms @ %.2f sleep} -> "
            "{%09.1f/%09.1fms, %03d iter, %.2e err, %03d early] %s\n",
            double(batchSize), double(affectedCount), sleepDuration, sleepProbability,
            ans.correctedTime, ans.time, ans.iterations, err, ear.iterations-1, technique
          );
        };
        // Find multi-threaded OpenMP-based Static PageRank (synchronous, no dead ends).
        auto a1 = pagerankBasicOmp(yt, init, {repeat}, fv);
        flog(a1, a1, "pagerankBasicOmp");
        // Find multi-threaded OpenMP-based Naive-dynamic PageRank (synchronous, no dead ends).
        auto a2 = pagerankBasicOmp(yt, &a0.ranks, {repeat}, fv);
        flog(a2, a1, "pagerankBasicNaiveDynamicOmp");
        // Find multi-threaded OpenMP-based Traversal-based Dynamic PageRank (synchronous, no dead ends).
        auto a3 = pagerankBasicDynamicTraversalOmp(x, xt, y, yt, deletions, insertions, &a0.ranks, {repeat}, fv);
        flog(a3, a1, "pagerankBasicDynamicTraversalOmp");
        // Find multi-threaded OpenMP-based Static Barrier-free PageRank (asynchronous, no dead ends).
        auto a4 = pagerankBarrierfreeOmp<true>(yt, init, {repeat}, fv);
        flog(a4, a1, "pagerankBarrierfreeOmp");
        // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free PageRank (asynchronous, no dead ends).
        auto a5 = pagerankBarrierfreeOmp<true>(yt, &a0.ranks, {repeat}, fv);
        flog(a5, a1, "pagerankBarrierfreeNaiveDynamicOmp");
        // Find multi-threaded OpenMP-based Traversal-based Dynamic Barrier-free PageRank (asynchronous, no dead ends).
        auto a6 = pagerankBarrierfreeDynamicTraversalOmp<true>(x, xt, y, yt, deletions, insertions, &a0.ranks, {repeat}, fv);
        flog(a6, a1, "pagerankBarrierfreeDynamicTraversalOmp");
      }
    }
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  OutDiGraph<uint32_t> x;
  readMtxOmpW(x, file); LOG(""); println(x);
  auto fl = [](auto u) { return true; };
  x = selfLoopOmp(x, None(), fl);      LOG(""); print(x);  printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegreeOmp(x); LOG(""); print(xt); printf(" (transposeWithDegree)\n");
  runExperiment(x, xt, repeat);
  printf("\n");
  return 0;
}
