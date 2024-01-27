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
#define MAX_THREADS 64
#endif
#ifndef REPEAT_BATCH
#define REPEAT_BATCH 1
#endif
#ifndef REPEAT_METHOD
#define REPEAT_METHOD 1
#endif




// PERFORM EXPERIMENT
// ------------------

template <class F>
inline void runThreads(F fn) {
  for (int t=NUM_THREADS_BEGIN; t<=NUM_THREADS_END; t NUM_THREADS_STEP) {
    omp_set_num_threads(t);
    fn(t);
    omp_set_num_threads(MAX_THREADS);
  }
}


template <class G, class F>
inline void runSleepFailures(const G& x, F fn) {
  // Randomly sleep after processing each vertex.
  for      (int d=FAILURE_DURATION_BEGIN;    d<=FAILURE_DURATION_END;    d FAILURE_DURATION_STEP) {
    for (double p=FAILURE_PROBABILITY_BEGIN; p<=FAILURE_PROBABILITY_END; p FAILURE_PROBABILITY_STEP) {
      for  (int t=FAILURE_THREADS_BEGIN;     t<=FAILURE_THREADS_END;     t FAILURE_THREADS_STEP) {
        chrono::milliseconds sd(d);
        double sp = p;
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
      double cp = p;
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
void runExperiment(G& x, H& xt, istream& fstream, size_t rows, size_t size, double batchFraction, size_t batchLength) {
  using  K = typename G::key_type;
  using  V = TYPE;
  vector<V> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat     = REPEAT_METHOD;
  int numThreads = MAX_THREADS;
  // Get ranks of vertices on original graph (static).
  auto fnop = [&](ThreadInfo *thread, auto v) {};
  auto r0   = pagerankBasicOmp(xt, init, {1, LI_NORM, 1e-100}, fnop);
  auto R10  = r0.ranks;
  auto R11  = r0.ranks;
  auto R20  = r0.ranks;
  auto R21  = r0.ranks;
  vector<tuple<K, K>> deletions;
  vector<tuple<K, K>> insertions;
  // Get ranks of vertices on updated graph (dynamic).
  for (int batchIndex=0; batchIndex<batchLength; ++batchIndex) {
    auto y = duplicate(x);
    insertions.clear();
    auto fb = [&](auto u, auto v, auto w) {
      insertions.push_back({u, v});
      y.addEdge(u, v);
    };
    readTemporalDo(fstream, false, false, rows, size_t(batchFraction * size), fb);
    updateOmpU(y);
    auto yt = transposeWithDegreeOmp(y);
    LOG(""); print(y); printf(" (insertions=%zu)\n", insertions.size());
    runThreads([&](int numThreads) {
      runFailures(y, [&](int failureDuration, double failureProbability, int failureThreads, auto fv) {
        // Follow a specific result logging format, which can be easily parsed later.
        auto flog  = [&](const auto& ans, const auto& ref, const char *technique) {
          auto err = liNormOmp(ans.ranks, ref.ranks);
          printf(
            "{-%.3e/+%.3e batchf, %04d batchi, %03d/%03d threads %04dms @ %.2e %s failure} -> "
            "{%09.1f/%09.1fms, %03d iter, %.2e err, %03d crashed] %s\n",
            0.0, batchFraction, batchIndex,
            failureThreads, numThreads, failureDuration, failureProbability, FAILURE_TYPE,
            ans.correctedTime, ans.time, ans.iterations, err, ans.crashedCount, technique
          );
        };
        auto s0 = pagerankBasicOmp(yt, init, {1, LI_NORM, 1e-100}, fnop);
        // Find multi-threaded OpenMP-based Static PageRank (synchronous, no dead ends).
        auto a0 = pagerankBasicOmp(yt, init, {repeat}, fv);
        flog(a0, s0, "pagerankBasicOmp");
        // Find multi-threaded OpenMP-based Naive-dynamic PageRank (synchronous, no dead ends).
        auto a1 = pagerankBasicOmp(yt, &R10, {repeat}, fv);
        flog(a1, s0, "pagerankBasicNaiveDynamicOmp");
        // Find multi-threaded OpenMP-based Frontier-based Dynamic PageRank (synchronous, no dead ends).
        auto a2 = pagerankBasicDynamicFrontierOmp(x, xt, y, yt, deletions, insertions, &R20, {repeat}, fv);
        flog(a2, s0, "pagerankBasicDynamicFrontierOmp");
        // Find multi-threaded OpenMP-based Static Barrier-free PageRank (asynchronous, no dead ends).
        auto b0 = pagerankBarrierfreeOmp<true>(yt, init, {repeat}, fv);
        flog(b0, s0, "pagerankBarrierfreeOmp");
        // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free PageRank (asynchronous, no dead ends).
        auto b1 = pagerankBarrierfreeOmp<true>(yt, &R11, {repeat}, fv);
        flog(b1, s0, "pagerankBarrierfreeNaiveDynamicOmp");
        // Find multi-threaded OpenMP-based Frontier-based Dynamic Barrier-free PageRank (asynchronous, no dead ends).
        auto b2 = pagerankBarrierfreeDynamicFrontierOmp<true>(x, xt, y, yt, deletions, insertions, &R21, {repeat}, fv);
        flog(b2, s0, "pagerankBarrierfreeDynamicFrontierOmp");
        // Update ranks.
        copyValuesOmpW(R10, a1.ranks);
        copyValuesOmpW(R20, a2.ranks);
        copyValuesOmpW(R11, b1.ranks);
        copyValuesOmpW(R21, b2.ranks);
      });
    });
    swap(x, y);
    swap(xt, yt);
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  size_t rows = strtoull(argv[2], nullptr, 10);
  size_t size = strtoull(argv[3], nullptr, 10);
  double batchFraction = strtod(argv[5], nullptr);
  size_t batchLength = strtoull(argv[6], nullptr, 10);
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  OutDiGraph<uint32_t> x;
  ifstream fstream(file);
  readTemporalOmpW(x, fstream, false, false, rows, size_t(0.90 * size)); LOG(""); print(x); printf(" (90%%)\n");
  auto fl = [](auto u) { return true; };
  x = selfLoopOmp(x, None(), fl);  LOG(""); print(x);  printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegreeOmp(x); LOG(""); print(xt); printf(" (transposeWithDegree)\n");
  runExperiment(x, xt, fstream, rows, size, batchFraction, batchLength);
  printf("\n");
  return 0;
}
