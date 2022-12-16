#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "pagerank.hxx"

using std::vector;
using std::swap;
using std::abs;




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
void pagerankFactor(vector<V>& a, const vector<K>& vdeg, V P, K i, K n) {
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
V pagerankTeleport(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  V a = (1-P)/N;
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}




// PAGERANK-CALCULATE
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
  a[v] = C0 + sumValuesAt(c, sliceIterable(xe, xv[v], xv[v+1]));
}


/**
 * Calculate rank for a given vertex, and get the change in rank value.
 * @param a current rank of each vertex (output)
 * @param r previous rank of each vertex
 * @param f rank scaling factor for each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param v given vertex
 * @param C0 common teleport rank contribution to each vertex
 * @returns change between previous and current rank value
 */
template <class K, class V>
inline V pagerankCalculateRankDelta(vector<V>& a, const vector<V>& r, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, K v, V C0) {
  V av = C0, rv = r[v];
  for (size_t i=xv[v], I=xv[v+1]; i<I; ++i) {
    K u = xe[i];
    av += r[u] * f[u];
  }
  a[v] = av;
  return abs(av - rv);
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
 */
template <class K, class V, class FV>
void pagerankCalculateRanks(vector<V>& a, const vector<V>& c, const vector<size_t>& xv, const vector<K>& xe, V C0, K i, K n, ThreadInfo *thread, FV fv) {
  for (K v=i; v<i+n; ++v)
  { pagerankCalculateRank(a, c, xv, xe, v, C0); fv(thread, v); }
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
inline V pagerankError(const vector<V>& x, const vector<V>& y, int EF, K i, K n) {
  switch (EF) {
    case 1:  return l1Norm(x, y, i, n);
    case 2:  return l2Norm(x, y, i, n);
    default: return liNorm(x, y, i, n);
  }
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
 * @returns pagerank result
 */
template <bool ASYNC=false, class H, class V, class KS, class NS, class FL, class FV>
PagerankResult<V> pagerankSeq(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const KS& ks, size_t i, const NS& ns, FL fl, FV fv) {
  using K  = typename H::key_type;
  K   N  = xt.order();
  V   P  = o.damping;
  V   E  = o.tolerance;
  int L  = o.maxIterations, l = 0;
  int EF = o.toleranceNorm;
  auto xv   = sourceOffsetsAs(xt, ks, size_t());
  auto xe   = destinationIndicesAs(xt, ks, K());
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
    pagerankFactor(f, vdeg, P, K(), N); multiplyValuesW(c, r, f, 0, N);  // calculate factors (f) and contributions (c)
    l = fl(e, ASYNC? r : a, r, c, f, xv, xe, vdeg, N, P, E, L, EF, K(i), ns, threads, fv);  // calculate ranks of vertices
    tcorrected += threadInfosMinDurationMilliseconds(threads, start);
  }, o.repeat);
  threadInfosDelete(threads);
  return {decompressContainer(xt, r, ks), l, t, tcorrected>0? tcorrected / o.repeat : t};
}
