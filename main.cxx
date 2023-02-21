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
void runExperiment(const G& x, const H& xt) {
  using  K = typename G::key_type;
  using  V = TYPE;
  vector<V> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat = REPEAT_METHOD;
  // Get ranks of vertices on original graph (static).
  auto fnop = [&](ThreadInfo *thread, auto v) {};
  auto a0   = pagerankBasicOmp(xt, init, {1}, fnop);
  auto b0   = pagerankBarrierfreeOmp<true>(xt, init, {1}, fnop);
  // Get ranks of vertices on updated graph (dynamic).
  runBatches(x, rnd, [&](const auto& y, const auto& yt, const auto& deletions, const auto& insertions) {
    vector<bool> vaff(max(x.span(), y.span()));
    pagerankAffectedTraversalW(vaff, x, y, deletions, insertions);
    size_t affectedCount = countIf(vaff, [](bool v) { return v==true; });
    runFailures(y, [&](int failureDuration, double failureProbability, int failureThreads, auto fv) {
      // Follow a specific result logging format, which can be easily parsed later.
      auto flog  = [&](const auto& ans, const auto& ref, const char *technique, int tfactor=1) {
        auto ear = pagerankBasicOmp(yt, &ans.ranks, {1}, fnop);
        auto err = l1NormOmp(ans.ranks, ref.ranks);
        LOG(
          "{-%.3e/+%.3e batch, %.3e aff, %03d/%03d threads %04dms @ %.2e %s failure} -> "
          "{%09.1f/%09.1fms, %03d iter, %.2e err, %03d early, %03d crashed] %s<tolerance_factor=%03d>\n",
          double(deletions.size()), double(insertions.size()), double(affectedCount),
          failureThreads, MAX_THREADS, failureDuration, failureProbability, FAILURE_TYPE,
          ans.correctedTime, ans.time, ans.iterations, err, ear.iterations-1, ans.crashedCount, technique, tfactor
        );
      };
      auto r1 = pagerankBasicOmp(yt, init, {1}, fnop);
      // Find multi-threaded OpenMP-based Static PageRank (synchronous, no dead ends).
      // auto a1 = pagerankBasicOmp(yt, init, {repeat}, fv);
      // flog(a1, r1, "pagerankBasicOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic PageRank (synchronous, no dead ends).
      auto a2 = pagerankBasicOmp(yt, &a0.ranks, {repeat}, fv);
      flog(a2, r1, "pagerankBasicNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Frontier-based Dynamic PageRank (synchronous, no dead ends).
      for (int T=1; T<256; T*=2) {
        auto a4 = pagerankBasicDynamicFrontierOmp<false, false, false>(x, xt, y, yt, deletions, insertions, &a0.ranks, {repeat, 1e-10/T}, fv);
        flog(a4, r1, "pagerankBasicDynamicFrontierOmp", T);
        auto a5 = pagerankBasicDynamicFrontierOmp<false, false, true> (x, xt, y, yt, deletions, insertions, &a0.ranks, {repeat, 1e-10/T}, fv);
        flog(a5, r1, "pagerankBasicDynamicFrontierOmpScaled", T);
      }
      // Find multi-threaded OpenMP-based Static Barrier-free PageRank (asynchronous, no dead ends).
      // auto b1 = pagerankBarrierfreeOmp<true>(yt, init, {repeat}, fv);
      // flog(b1, r1, "pagerankBarrierfreeOmp");
      // Find multi-threaded OpenMP-based Naive-dynamic Barrier-free PageRank (asynchronous, no dead ends).
      auto b2 = pagerankBarrierfreeOmp<true>(yt, &b0.ranks, {repeat}, fv);
      flog(b2, r1, "pagerankBarrierfreeNaiveDynamicOmp");
      // Find multi-threaded OpenMP-based Frontier-based Dynamic Barrier-free PageRank (asynchronous, no dead ends).
      for (int T=1; T<256; T*=2) {
        auto b4 = pagerankBarrierfreeDynamicFrontierOmp<true, false, false> (x, xt, y, yt, deletions, insertions, &b0.ranks, {repeat, 1e-10/T}, fv);
        flog(b4, r1, "pagerankBarrierfreeDynamicFrontierOmp", T);
        auto b5 = pagerankBarrierfreeDynamicFrontierOmp<true, false, true>  (x, xt, y, yt, deletions, insertions, &b0.ranks, {repeat, 1e-10/T}, fv);
        flog(b5, r1, "pagerankBarrierfreeDynamicFrontierOmpScaled", T);
      }
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
