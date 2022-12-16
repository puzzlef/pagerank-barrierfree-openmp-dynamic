#pragma once
#include <utility>
#include <chrono>
#include <random>
#include <atomic>
#include <vector>
#include "_main.hxx"

using std::chrono::system_clock;
using std::random_device;
using std::default_random_engine;
using std::atomic;
using std::vector;
using std::move;




// PAGERANK-OPTIONS
// ----------------

enum NormFunction {
  L0_NORM = 0,
  L1_NORM = 1,
  L2_NORM = 2,
  LI_NORM = 3
};


template <class V>
struct PagerankOptions {
  int repeat;
  int toleranceNorm;
  V   tolerance;
  V   damping;
  int maxIterations;

  PagerankOptions(int repeat=1, int toleranceNorm=LI_NORM, V tolerance=1e-10, V damping=0.85, int maxIterations=500) :
  repeat(repeat), toleranceNorm(toleranceNorm), tolerance(tolerance), damping(damping), maxIterations(maxIterations) {}
};




// PAGERANK-RESULT
// ---------------

template <class V>
struct PagerankResult {
  vector<V> ranks;
  int   iterations;
  float time;
  float correctedTime;

  PagerankResult() :
  ranks(), iterations(0), time(0), correctedTime(0) {}

  PagerankResult(vector<V>&& ranks, int iterations=0, float time=0, float correctedTime=0) :
  ranks(ranks), iterations(iterations), time(time), correctedTime(correctedTime) {}

  PagerankResult(vector<V>& ranks, int iterations=0, float time=0, float correctedTime=0) :
  ranks(move(ranks)), iterations(iterations), time(time), correctedTime(correctedTime) {}
};




// THREAD-INFO
// -----------

struct ThreadInfo {
  random_device dev;          // used for random sleeps
  default_random_engine rnd;  // used for random sleeps
  system_clock::time_point stop;  // stop time point
  int id;                     // thread number
  int iteration;              // current iteration

  ThreadInfo(int id) :
  dev(), rnd(dev()), id(id), iteration(0), stop() {}
  inline void clear() { iteration = 0; stop = system_clock::time_point(); }
};


inline vector<ThreadInfo*> threadInfos(int N) {
  vector<ThreadInfo*> threads(N);
  for (int i=0; i<N; ++i)
    threads[i] = new ThreadInfo(i);
  return threads;
}

inline void threadInfosDelete(const vector<ThreadInfo*>& threads) {
  int N = threads.size();
  for (int i=0; i<N; ++i)
    delete threads[i];
}

inline void threadInfosClear(const vector<ThreadInfo*>& threads) {
  int N = threads.size();
  for (int i=0; i<N; ++i)
    threads[i]->clear();
}

inline float threadInfosMinDurationMilliseconds(const vector<ThreadInfo*>& threads, system_clock::time_point start) {
  int N = threads.size(); float a = 0;
  for (int i=0; i<N; ++i) {
    if (threads[i]->stop <= start) continue;
    float t = durationMilliseconds(start, threads[i]->stop);
    if (a==0 || a>t) a = t;
  }
  return a;
}
