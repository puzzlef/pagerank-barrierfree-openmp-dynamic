#pragma once
#include <algorithm>
#include <vector>
#include <omp.h>
#include "_main.hxx"
#include "transpose.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"

using std::vector;
using std::swap;




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
V pagerankBarrierfreeTeleportOmp(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  V a = (1-P)/N;
  #pragma omp for schedule(auto) reduction(+:a) nowait
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}




// PAGERANK-CALCULATE
// ------------------
// For rank calculation from in-edges.

/**
 * Calculate ranks for vertices in a graph.
 * @param e change in rank for each vertex below tolerance? (updated)
 * @param a current rank of each vertex (updated)
 * @param r previous rank of each vertex
 * @param f rank scaling factor for each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param C0 common teleport rank contribution to each vertex
 * @param E tolerance [10^-10]
 * @param i vertex start
 * @param n vertex count
 * @param thread information on current thread (updated)
 * @param fv per vertex processing (thread, vertex)
 */
template <class K, class V, class FV>
void pagerankBarrierfreeCalculateRanksOmp(vector<int>& e, vector<V>& a, const vector<V>& r, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, V C0, V E, K i, K n, ThreadInfo *thread, FV fv) {
  #pragma omp for schedule(dynamic, 2048) nowait
  for (K v=i; v<i+n; ++v) {
    V ev = pagerankCalculateRankDelta(a, r, f, xv, xe, v, C0);
    if (ev<=E && e[v]==0) e[v] = 1;  // LI_NORM
    fv(thread, v);
  }
}




// PAGERANK-CONVERGED
// ------------------
// For convergence check.

/**
 * Check if ranks of all vertices have converged.
 * @param e change in rank for each vertex below tolerance?
 * @param i vertex start
 * @param n vertex count
 * @returns ranks converged?
 */
template <class K>
inline bool pagerankBarrierfreeConverged(const vector<int>& e, K i, K n) {
  for (K u=i; u<i+n; ++u)
    if (!e[u]) return false;
  return true;
}




// PAGERANK-LOOP
// -------------

template <bool ASYNC=false, bool DEAD=false, class K, class V, class FV>
int pagerankBarrierfreeOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, vector<V>& c, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, const vector<K>& vdeg, K N, V P, V E, int L, int EF, K i, K n, vector<ThreadInfo*>& threads, FV fv) {
  if (EF!=LI_NORM) return 0;
  #pragma omp parallel
  {
    int  t = omp_get_thread_num();
    int& l = threads[t]->iteration;
    while (l<L) {
      V C0 = DEAD? pagerankBarrierfreeTeleportOmp(r, vdeg, P, N) : (1-P)/N;
      pagerankBarrierfreeCalculateRanksOmp(e, a, r, f, xv, xe, C0, E, i, n, threads[t], fv); ++l;  // update ranks of vertices
      if (!ASYNC) swap(a, r);                            // final ranks in (r)
      if (pagerankBarrierfreeConverged(e, i, n)) break;  // check tolerance
    }
    threads[t]->stop = timeNow();
  }
  int l = 0;
  for (int t=0; t<threads.size(); ++t)
    if (threads[t]->iteration > l) l = threads[t]->iteration;
  if (!ASYNC && (l & 1)==1) swap(a, r);
  return l;
}




// STATIC/NAIVE-DYNAMIC PAGERANK
// -----------------------------

/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @param fv per vertex processing (thread, vertex)
 * @returns pagerank result
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V, class FV>
inline PagerankResult<V> pagerankBarrierfreeOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  using K = typename H::key_type;
  K    N  = xt.order();  if (N==0) return {};
  auto ks = vertexKeys(xt);
  return pagerankOmp<ASYNC>(xt, q, o, ks, 0, N, pagerankBarrierfreeOmpLoop<ASYNC, DEAD, K, V, FV>, fv);
}

/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @returns pagerank result
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V>
inline PagerankResult<V> pagerankBarrierfreeOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o) {
  auto fv = [](ThreadInfo *thread, auto v) {};
  return pagerankBarrierfreeOmp<ASYNC, DEAD>(xt, q, o, fv);
}
