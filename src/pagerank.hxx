#pragma once
#include <utility>
#include <chrono>
#include <random>
#include <atomic>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "csr.hxx"
#include "vertices.hxx"
#include "transpose.hxx"
#include "dfs.hxx"
#include "components.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::random_device;
using std::default_random_engine;
using std::chrono::system_clock;
using std::tuple;
using std::vector;
using std::atomic;
using std::move;
using std::max;




// PAGERANK OPTIONS
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




// PAGERANK RESULT
// ---------------

template <class V>
struct PagerankResult {
  vector<V> ranks;
  int   iterations;
  float time;
  float correctedTime;
  int   crashedCount;

  PagerankResult() :
  ranks(), iterations(0), time(0), correctedTime(0), crashedCount(0) {}

  PagerankResult(vector<V>&& ranks, int iterations=0, float time=0, float correctedTime=0, int crashedCount=0) :
  ranks(ranks), iterations(iterations), time(time), correctedTime(correctedTime), crashedCount(crashedCount) {}

  PagerankResult(vector<V>& ranks, int iterations=0, float time=0, float correctedTime=0, int crashedCount=0) :
  ranks(move(ranks)), iterations(iterations), time(time), correctedTime(correctedTime), crashedCount(crashedCount) {}
};




// PAGERANK DATA
// -------------
// Using Pagerank Data for performance!

template <class G>
struct PagerankData {
  using K = typename G::key_type;
  vector2d<K> components;
  G blockgraph;
  G blockgraphTranspose;
};

template <class G, class K>
auto blockgraphD(const G& x, const vector2d<K>& cs, const PagerankData<G> *D) {
  return D? D->blockgraph : blockgraph(x, cs);
}

template <class G>
auto blockgraphTransposeD(const G& b, const PagerankData<G> *D) {
  return D? D->blockgraphTranspose : transpose(b);
}

template <class G, class H>
auto componentsD(const G& x, const H& xt, const PagerankData<G> *D) {
  return D? D->components : components(x, xt);
}




// THREAD INFO
// -----------

struct ThreadInfo {
  random_device dev;              // used for random sleeps
  default_random_engine rnd;      // used for random sleeps
  system_clock::time_point stop;  // stop time point
  int   id;         // thread number
  float iteration;  // current iteration
  bool  crashed;    // error occurred?

  ThreadInfo(int id) :
  dev(), rnd(dev()), stop(), id(id), iteration(0), crashed(false) {}
  inline void clear() { stop = system_clock::time_point(); iteration = 0; crashed = false; }
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

inline int threadInfosCrashedCount(const vector<ThreadInfo*>& threads) {
  int N = threads.size(), a = 0;
  for (int i=0; i<N; ++i)
    if (threads[i]->crashed) ++a;
  return a;
}

inline float threadInfosMaxIteration(const vector<ThreadInfo*>& threads) {
  int N = threads.size(); float a = 0;
  for (int i=0; i<N; ++i) {
    if (threads[i]->crashed) continue;
    if (threads[i]->iteration > a) a = threads[i]->iteration;
  }
  return a;
}

inline float threadInfosMinDuration(const vector<ThreadInfo*>& threads, system_clock::time_point start) {
  int N = threads.size(); float a = 0;
  for (int i=0; i<N; ++i) {
    if (threads[i]->stop <= start || threads[i]->crashed) continue;
    float t = duration(start, threads[i]->stop);
    if (a==0 || a>t) a = t;
  }
  return a;
}




// PAGERANK FACTOR
// ---------------
// For contribution factors of vertices (unchanging).

/**
 * Calculate rank scaling factor for each vertex.
 * @param a rank scaling factor for each vertex (output)
 * @param vdeg out-degree of each vertex
 * @param P damping factor [0.85]
 * @param i vertex start
 * @param n vertex count
 */
template <class K, class V>
inline void pagerankFactor(vector<V>& a, const vector<K>& vdeg, V P, K i, K n) {
  for (K u=i; u<i+n; ++u) {
    K  d = vdeg[u];
    a[u] = d>0? P/d : 0;
  }
}


#ifdef OPENMP
template <class K, class V>
inline void pagerankFactorOmp(vector<V>& a, const vector<K>& vdeg, V P, K i, K n) {
  #pragma omp parallel for schedule(auto)
  for (K u=i; u<i+n; ++u) {
    K  d = vdeg[u];
    a[u] = d>0? P/d : 0;
  }
}
#endif




// PAGERANK TELEPORT
// -----------------
// For teleport contribution from vertices (inc. dead ends).

/**
 * Find total teleport contribution from each vertex (inc. deade ends).
 * @param r rank of each vertex
 * @param vdeg out-degree of each vertex
 * @param P damping factor [0.85]
 * @param N total number of vertices
 * @returns common teleport rank contribution to each vertex
 */
template <class K, class V>
inline V pagerankTeleport(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  V a = (1-P)/N;
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}


#ifdef OPENMP
template <class K, class V>
inline V pagerankTeleportOmp(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  V a = (1-P)/N;
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}
#endif




// PAGERANK CALCULATE
// ------------------
// For rank calculation from in-edges.

/**
 * Calculate rank for a given vertex.
 * @param a current rank of each vertex (output)
 * @param c rank contribution from each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param v given vertex
 * @param C0 common teleport rank contribution to each vertex
 */
template <class K, class V>
inline void pagerankCalculateRank(vector<V>& a, const vector<V>& c, const vector<size_t>& xv, const vector<K>& xe, K v, V C0) {
  V av = C0;
  for (size_t i=xv[v], I=xv[v+1]; i<I; ++i)
    av += c[xe[i]];
  a[v] = av;
}


/**
 * Calculate ranks for vertices in a graph.
 * @param a current rank of each vertex (output)
 * @param c rank contribution from each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param C0 common teleport rank contribution to each vertex
 * @param i vertex start
 * @param n vertex count
 * @param thread information on current thread (updated)
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 */
template <class K, class V, class FV, class FA>
inline void pagerankCalculateRanks(vector<V>& a, const vector<V>& c, const vector<size_t>& xv, const vector<K>& xe, V C0, K i, K n, ThreadInfo *thread, FV fv, FA fa) {
  for (K v=i; v<i+n; ++v) {
    if (!fa(v)) continue;
    pagerankCalculateRank(a, c, xv, xe, v, C0);
    fv(thread, v);
  }
}


#ifdef OPENMP
template <class K, class V, class FV, class FA>
inline void pagerankCalculateRanksOmp(vector<V>& a, const vector<V>& c, const vector<size_t>& xv, const vector<K>& xe, V C0, K i, K n, vector<ThreadInfo*>& threads, FV fv, FA fa) {
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K v=i; v<i+n; ++v) {
    if (!fa(v)) continue;
    int t = omp_get_thread_num();
    pagerankCalculateRank(a, c, xv, xe, v, C0);
    fv(threads[t], v);
  }
}
#endif




// PAGERANK ERROR
// --------------
// For convergence check.

/**
 * Get the error between two rank vectors.
 * @param x first rank vector
 * @param y second rank vector
 * @param EF error function (L1/L2/LI)
 * @param i vertex start
 * @param n vertex count
 * @returns error between the two rank vectors
 */
template <class K, class V>
inline V pagerankError(const vector<V>& x, const vector<V>& y, int EF, K i, K n) {
  switch (EF) {
    case 1:  return l1Norm(x, y, i, n);
    case 2:  return l2Norm(x, y, i, n);
    default: return liNorm(x, y, i, n);
  }
}


#ifdef OPENMP
template <class K, class V>
inline V pagerankErrorOmp(const vector<V>& x, const vector<V>& y, int EF, K i, K N) {
  switch (EF) {
    case 1:  return l1NormOmp(x, y, i, N);
    case 2:  return l2NormOmp(x, y, i, N);
    default: return liNormOmp(x, y, i, N);
  }
}
#endif




// PAGERANK AFFECTED (TRAVERSAL)
// -----------------------------

/**
 * Find affected vertices due to a batch update.
 * @param x original graph
 * @param y updated graph
 * @param ft is vertex affected? (u)
 * @returns affected flags
 */
template <class G, class FT>
inline auto pagerankAffectedTraversal(const G& x, const G& y, FT ft) {
  auto fn = [](auto u) {};
  vector<bool> vis(max(x.span(), y.span()));
  y.forEachVertexKey([&](auto u) {
    if (!ft(u)) return;
    dfsVisitedForEachW(vis, x, u, fn);
    dfsVisitedForEachW(vis, y, u, fn);
  });
  return vis;
}


/**
 * Find affected vertices due to a batch update.
 * @param y original graph
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @returns affected flags
 */
template <class G, class K>
inline auto pagerankAffectedTraversal(const G& x, const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions) {
  auto fn = [](K u) {};
  vector<bool> vis(max(x.span(), y.span()));
  for (const auto& [u, v] : deletions)
    dfsVisitedForEachW(vis, x, u, fn);
  for (const auto& [u, v] : insertions)
    dfsVisitedForEachW(vis, y, u, fn);
  return vis;
}




// PAGERANK-SEQ
// ------------
// For single-threaded (sequential) PageRank implementation.

/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @param ks vertices (keys) to process
 * @param i vertex start
 * @param ns vertex count(s)
 * @param fl update loop
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 * @returns pagerank result
 */
template <bool ASYNC=false, class H, class V, class KS, class NS, class FL, class FV, class FA>
PagerankResult<V> pagerankSeq(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const KS& ks, size_t i, const NS& ns, FL fl, FV fv, FA fa) {
  using K  = typename H::key_type;
  K   N  = xt.order();
  V   P  = o.damping;
  V   E  = o.tolerance;
  int L  = o.maxIterations, l = 0;
  int EF = o.toleranceNorm;
  auto xv   = sourceOffsets(xt, ks);
  auto xe   = destinationIndices(xt, ks);
  auto vdeg = vertexData(xt, ks);
  vector<ThreadInfo*> threads = threadInfos(1);
  vector<int> e(N); vector<V> a(N), r(N), c(N), f(N), qc;
  if (q) qc = compressContainer(xt, *q, ks);
  float tcorrected = 0;
  float t = measureDuration([&]() {
    auto start = timeNow();
    threadInfosClear(threads);
    fillValueU(e, 0);
    if (q) copyValuesW(r, qc);
    else   fillValueU (r, V(1)/N);
    if (!ASYNC) copyValuesW(a, r);
    pagerankFactor(f, vdeg, P, K(), N); multiplyValuesW(c, r, f, 0, N);  // calculate factors (f) and contributions (c)
    l = fl(e, ASYNC? r : a, r, c, f, xv, xe, vdeg, N, P, E, L, EF, K(i), ns, threads, fv, fa);  // calculate ranks of vertices
    tcorrected += threadInfosMinDuration(threads, start);
  }, o.repeat);
  float correctedTime = tcorrected>0? tcorrected / o.repeat : t;
  int   crashedCount  = threadInfosCrashedCount(threads);
  threadInfosDelete(threads);
  return {decompressContainer(xt, r, ks), l, t, correctedTime, crashedCount};
}




// PAGERANK-OMP
// ------------
// For multi-threaded OpenMP-based PageRank implementation.

#ifdef OPENMP
/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @param ks vertices (keys) to process
 * @param i vertex start
 * @param ns vertex count(s)
 * @param fl update loop
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 * @returns pagerank result
 */
template <bool ASYNC=false, class H, class V, class KS, class NS, class FL, class FV, class FA>
PagerankResult<V> pagerankOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const KS& ks, size_t i, const NS& ns, FL fl, FV fv, FA fa) {
  using K  = typename H::key_type;
  K   N  = xt.order();
  V   P  = o.damping;
  V   E  = o.tolerance;
  int L  = o.maxIterations, l = 0;
  int EF = o.toleranceNorm;
  int TH = omp_get_max_threads();
  auto xv   = sourceOffsets(xt, ks);
  auto xe   = destinationIndices(xt, ks);
  auto vdeg = vertexData(xt, ks);
  vector<ThreadInfo*> threads = threadInfos(TH);
  vector<int> e(N); vector<V> a(N), r(N), c(N), f(N), qc;
  if (q) qc = compressContainer(xt, *q, ks);
  float tcorrected = 0;
  float t = measureDuration([&]() {
    auto start = timeNow();
    threadInfosClear(threads);
    fillValueU(e, 0);
    if (q) copyValuesOmpW(r, qc);
    else   fillValueOmpU (r, V(1)/N);
    if (!ASYNC) copyValuesOmpW(a, r);
    pagerankFactorOmp(f, vdeg, P, K(), N); multiplyValuesOmpW(c, r, f, 0, N);  // calculate factors (f) and contributions (c)
    l = fl(e, ASYNC? r : a, r, c, f, xv, xe, vdeg, N, P, E, L, EF, K(i), ns, threads, fv, fa);  // calculate ranks of vertices
    tcorrected += threadInfosMinDuration(threads, start);
  }, o.repeat);
  float correctedTime = tcorrected>0? tcorrected / o.repeat : t;
  int   crashedCount  = threadInfosCrashedCount(threads);
  threadInfosDelete(threads);
  return {decompressContainer(xt, r, ks), l, t, correctedTime, crashedCount};
}
#endif
