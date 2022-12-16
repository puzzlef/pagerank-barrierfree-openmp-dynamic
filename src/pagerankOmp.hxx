#pragma once
#include <algorithm>
#include <vector>
#include <omp.h>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"

using std::vector;
using std::swap;
using std::max;




// PAGERANK-FACTOR
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
void pagerankFactorOmp(vector<V>& a, const vector<K>& vdeg, V P, K i, K n) {
  if (n < SIZE_MIN_OMPM) { pagerankFactor(a, vdeg, P, i, n); return; }
  #pragma omp parallel for schedule(auto)
  for (K u=i; u<i+n; ++u) {
    K  d = vdeg[u];
    a[u] = d>0? P/d : 0;
  }
}




// PAGERANK-TELEPORT
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
V pagerankTeleportOmp(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  if (N < SIZE_MIN_OMPR) return pagerankTeleport(r, vdeg, P, N);
  V a = (1-P)/N;
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}




// PAGERANK-CALCULATE
// ------------------
// For rank calculation from in-edges.

/**
 * Calculate ranks for vertices in a graph.
 * @param a current rank of each vertex (output)
 * @param c rank contribution from each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param C0 common teleport rank contribution to each vertex
 * @param i vertex start
 * @param n vertex count
 * @param threads information on each thread (updated)
 * @param fv per vertex processing (thread, vertex)
 */
template <class K, class V, class FV>
void pagerankCalculateRanksOmp(vector<V>& a, const vector<V>& c, const vector<size_t>& xv, const vector<K>& xe, V C0, K i, K n, vector<ThreadInfo*>& threads, FV fv) {
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K v=i; v<i+n; ++v) {
    int t = omp_get_thread_num();
    pagerankCalculateRank(a, c, xv, xe, v, C0); fv(threads[t], v);
  }
}




// PAGERANK-ERROR
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
inline V pagerankErrorOmp(const vector<V>& x, const vector<V>& y, int EF, K i, K N) {
  switch (EF) {
    case 1:  return l1NormOmp(x, y, i, N);
    case 2:  return l2NormOmp(x, y, i, N);
    default: return liNormOmp(x, y, i, N);
  }
}




// PAGERANK-OMP
// ------------
// For multi-threaded OpenMP-based PageRank implementation.

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
 * @returns pagerank result
 */
template <bool ASYNC=false, class H, class V, class KS, class NS, class FL, class FV>
PagerankResult<V> pagerankOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const KS& ks, size_t i, const NS& ns, FL fl, FV fv) {
  using K  = typename H::key_type;
  K   N  = xt.order();
  V   P  = o.damping;
  V   E  = o.tolerance;
  int L  = o.maxIterations, l = 0;
  int EF = o.toleranceNorm;
  int TH = omp_get_max_threads();
  auto xv   = sourceOffsetsAs(xt, ks, size_t());
  auto xe   = destinationIndicesAs(xt, ks, K());
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
    pagerankFactorOmp(f, vdeg, P, K(), N); multiplyValuesOmpW(c, r, f, 0, N);  // calculate factors (f) and contributions (c)
    l = fl(e, ASYNC? r : a, r, c, f, xv, xe, vdeg, N, P, E, L, EF, K(i), ns, threads, fv);  // calculate ranks of vertices
    tcorrected += threadInfosMinDurationMilliseconds(threads, start);
  }, o.repeat);
  threadInfosDelete(threads);
  return {decompressContainer(xt, r, ks), l, t, tcorrected>0? tcorrected / o.repeat : t};
}
