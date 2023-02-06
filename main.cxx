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




// Fixed config
#ifndef TYPE
#define TYPE double
#endif
#ifndef MAX_THREADS
#define MAX_THREADS 32
#endif
#ifndef REPEAT_BATCH
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
#define REPEAT_METHOD 5
#endif




// GENERATE BATCH
// --------------

template <class G, class R>
inline auto addRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v);
    insertions.push_back(make_tuple(u, v));
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(a, rnd, i, n, None(), fe); }, retries);
  updateOmpU(a);
  return insertions;
}


template <class G, class R>
inline auto removeRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    deletions.push_back(make_tuple(u, v));
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(a, rnd, i, n, fe); }, retries);
  updateOmpU(a);
  return deletions;
}




// PERFORM EXPERIMENT
// ------------------

template <class G, class R, class F>
inline void runAbsoluteBatches(const G& x, R& rnd, F fn) {
  size_t d = BATCH_DELETIONS_BEGIN;
  size_t i = BATCH_INSERTIONS_BEGIN;
  while (true) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      auto deletions  = removeRandomEdges(y, rnd, d, 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, i, 1, x.span()-1);
      auto yt = transposeWithDegreeOmp(y);
      fn(y, yt, deletions, insertions);
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, size_t(BATCH_DELETIONS_END));
    i = min(i, size_t(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runRelativeBatches(const G& x, R& rnd, F fn) {
  double d = BATCH_DELETIONS_BEGIN;
  double i = BATCH_INSERTIONS_BEGIN;
  while (true) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      auto deletions  = removeRandomEdges(y, rnd, size_t(d * x.size()), 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, size_t(i * x.size()), 1, x.span()-1);
      auto yt = transposeWithDegreeOmp(y);
      fn(y, yt, deletions, insertions);
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, double(BATCH_DELETIONS_END));
    i = min(i, double(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runBatches(const G& x, R& rnd, F fn) {
  if (BATCH_UNIT=="%") runRelativeBatches(x, rnd, fn);
  else runAbsoluteBatches(x, rnd, fn);
}


template <class G, class F>
inline void runSleepFailures(const G& x, F fn) {
  // Randomly sleep after processing each vertex.
  for      (int d=FAILURE_DURATION_BEGIN;    d<=FAILURE_DURATION_END;    d FAILURE_DURATION_STEP) {
    for (double p=FAILURE_PROBABILITY_BEGIN; p<=FAILURE_PROBABILITY_END; p FAILURE_PROBABILITY_STEP) {
      for  (int t=FAILURE_THREADS_BEGIN;     t<=FAILURE_THREADS_END;     t FAILURE_THREADS_STEP) {
        chrono::milliseconds sd(d);
        double sp = p / x.order();
        auto  fv = [&](ThreadInfo *thread, auto v) {
          uniform_real_distribution<double> dis(0.0, 1.0);
          if (thread->id < t && dis(thread->rnd) < sp) this_thread::sleep_for(sd);
        };
        fn(d, p, t, fv);
      }
    }
  }
}


template <class G, class F>
inline void runCrashFailures(const G& x, F fn) {
  // Randomly crash (simulated) after processing each vertex.
  for (double p=FAILURE_PROBABILITY_BEGIN; p<=FAILURE_PROBABILITY_END; p FAILURE_PROBABILITY_STEP) {
    for  (int t=FAILURE_THREADS_BEGIN;     t<=FAILURE_THREADS_END;     t FAILURE_THREADS_STEP) {
      double cp = p / x.order();
      auto  fv = [&](ThreadInfo *thread, auto v) {
        uniform_real_distribution<double> dis(0.0, 1.0);
        if (thread->id < t && dis(thread->rnd) < cp) thread->crashed = true;
      };
      fn(0, p, t, fv);
    }
  }
}


template <class G, class F>
inline void runFailures(const G& x, F fn) {
  auto fnop = [](ThreadInfo *thread, auto v) {};
  if      (FAILURE_TYPE=="sleep") runSleepFailures(x, fn);
  else if (FAILURE_TYPE=="crash") runCrashFailures(x, fn);
  else fn(0, 0.0f, 0, fnop);
}


template <class G, class H>
void runExperiment(const G& x, const H& xt) {
  using  K = typename G::key_type;
  using  V = TYPE;
  vector<V> *init = nullptr;
  PagerankData<G> *nodata = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat = REPEAT_METHOD;
  // Get ranks of vertices on original graph (static).
  auto fnop = [&](ThreadInfo *thread, auto v) {};
  auto a0   = pagerankBasicOmp(xt, init, {1}, fnop);
  auto b0   = pagerankBarrierfreeOmp<true>(xt, init, {1}, fnop);
  auto c0   = pagerankBarrierfreeMonolithicOmp<true>(x, xt, init, {1}, nodata, fnop);
  auto d0   = pagerankBarrierfreeLevelwiseOmp<true> (x, xt, init, {1}, nodata, fnop);
  // Get ranks of vertices on updated graph (dynamic).
  runBatches(x, rnd, [&](const auto& y, const auto& yt, const auto& deletions, const auto& insertions) {
    // Generate common Pagerank data for speeding up experiment.
    auto cs = components(y, yt);
    auto b  = blockgraph(y, cs);
    auto bt = transpose(b);
    auto gs = levelwiseGroupedComponentsFrom(cs, b, bt);
    size_t CS = cs.size(), GS = gs.size();
    PagerankData<G> C {move(cs), move(b), move(bt)};
    size_t affectedCount = countValue(pagerankAffectedTraversal(x, y, deletions, insertions), true);
    runFailures(y, [&](int failureDuration, double failureProbability, int failureThreads, auto fv) {
      // Follow a specific result logging format, which can be easily parsed later.
      auto flog  = [&](const auto& ans, const auto& ref, const char *technique) {
        auto ear = pagerankBasicOmp(yt, &ans.ranks, {1}, fnop);
        auto err = l1NormOmp(ans.ranks, ref.ranks);
        LOG(
          "{-%.3e/+%.3e batch, %.3e aff, %.3e sccs, %.3e levels, %03d/%03d threads %04dms @ %.2e %s failure} -> "
          "{%09.1f/%09.1fms, %03d iter, %.2e err, %03d early] %s\n",
          double(deletions.size()), double(insertions.size()), double(affectedCount), double(CS), double(GS),
          failureThreads, MAX_THREADS, failureDuration, failureProbability, FAILURE_TYPE,
          ans.correctedTime, ans.time, ans.iterations, err, ear.iterations-1, technique
        );
      };
      auto r1 = pagerankBasicOmp(yt, init, {1}, fnop);
      // Find multi-threaded OpenMP-based Static PageRank (synchronous, no dead ends).
      auto a1 = pagerankBasicOmp(yt, init, {repeat}, fv);
      flog(a1, r1, "pagerankBasicOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic PageRank (synchronous, no dead ends).
      auto a2 = pagerankBasicOmp(yt, &a0.ranks, {repeat}, fv);
      flog(a2, r1, "pagerankBasicNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Traversal-based Dynamic PageRank (synchronous, no dead ends).
      auto a3 = pagerankBasicDynamicTraversalOmp(x, xt, y, yt, deletions, insertions, &a0.ranks, {repeat}, fv);
      flog(a3, r1, "pagerankBasicDynamicTraversalOmp");
      // Find multi-threaded OpenMP-based Static Barrier-free PageRank (asynchronous, no dead ends).
      auto b1 = pagerankBarrierfreeOmp<true>(yt, init, {repeat}, fv);
      flog(b1, r1, "pagerankBarrierfreeOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free PageRank (asynchronous, no dead ends).
      auto b2 = pagerankBarrierfreeOmp<true>(yt, &b0.ranks, {repeat}, fv);
      flog(b2, r1, "pagerankBarrierfreeNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Traversal-based Dynamic Barrier-free PageRank (asynchronous, no dead ends).
      auto b3 = pagerankBarrierfreeDynamicTraversalOmp<true>(x, xt, y, yt, deletions, insertions, &b0.ranks, {repeat}, fv);
      flog(b3, r1, "pagerankBarrierfreeDynamicTraversalOmp");
      // Find multi-threaded OpenMP-based Static Barrier-free Monolithic PageRank (asynchronous, no dead ends).
      auto c1 = pagerankBarrierfreeMonolithicOmp<true>(y, yt, init, {repeat}, &C, fv);
      flog(c1, r1, "pagerankBarrierfreeMonolithicOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free Monolithic PageRank (asynchronous, no dead ends).
      auto c2 = pagerankBarrierfreeMonolithicOmp<true>(y, yt, &c0.ranks, {repeat}, &C, fv);
      flog(c2, r1, "pagerankBarrierfreeMonolithicNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Traversal-based Dynamic Barrier-free Monolithic PageRank (asynchronous, no dead ends).
      auto c3 = pagerankBarrierfreeMonolithicDynamicTraversalOmp<true>(x, xt, y, yt, deletions, insertions, &c0.ranks, {repeat}, &C, fv);
      flog(c3, r1, "pagerankBarrierfreeMonolithicDynamicTraversalOmp");
      // Find multi-threaded OpenMP-based Static Barrier-free Levelwise PageRank (asynchronous, no dead ends).
      auto d1 = pagerankBarrierfreeLevelwiseOmp<true>(y, yt, init, {repeat}, &C, fv);
      flog(d1, r1, "pagerankBarrierfreeLevelwiseOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free Levelwise PageRank (asynchronous, no dead ends).
      auto d2 = pagerankBarrierfreeLevelwiseOmp<true>(y, yt, &d0.ranks, {repeat}, &C, fv);
      flog(d2, r1, "pagerankBarrierfreeLevelwiseNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Traversal-based Dynamic Barrier-free Levelwise PageRank (asynchronous, no dead ends).
      auto d3 = pagerankBarrierfreeLevelwiseDynamicTraversalOmp<true>(x, xt, y, yt, deletions, insertions, &d0.ranks, {repeat}, &C, fv);
      flog(d3, r1, "pagerankBarrierfreeLevelwiseDynamicTraversalOmp");
    });
  });
}


int main(int argc, char **argv) {
  char *file = argv[1];
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  OutDiGraph<uint32_t> x;
  readMtxOmpW(x, file); LOG(""); println(x);
  auto fl = [](auto u) { return true; };
  x = selfLoopOmp(x, None(), fl);      LOG(""); print(x);  printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegreeOmp(x); LOG(""); print(xt); printf(" (transposeWithDegree)\n");
  runExperiment(x, xt);
  printf("\n");
  return 0;
}
