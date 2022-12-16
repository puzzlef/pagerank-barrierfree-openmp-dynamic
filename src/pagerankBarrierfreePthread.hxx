#pragma once
#include <algorithm>
#include <atomic>
#include <vector>
#include <pthread.h>
#include "_main.hxx"
#include "transpose.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"
#include "pagerankBarrierfreeOmp.hxx"

using std::atomic;
using std::vector;
using std::swap;
using std::min;




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
inline void pagerankBarrierfreeCalculateRanksPthread(vector<int>& e, vector<V>& a, const vector<V>& r, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, V C0, V E, K i, K n, ThreadInfo *thread, FV fv) {
  const K CHUNK_SIZE = 2048;
  parallelForDynamicNoWaitPthread(i, i+n, CHUNK_SIZE, [&](K v) {
    V ev = pagerankCalculateRankDelta(a, r, f, xv, xe, v, C0);
    if (ev<=E && e[v]==0) e[v] = 1;  // LI_NORM
    fv(thread, v);
  });
}




// PAGERANK-LOOP
// -------------

template <bool ASYNC=false, bool DEAD=false, class K, class V, class FV>
inline int pagerankBarrierfreePthreadLoop(vector<int>& e, vector<V>& a, vector<V>& r, vector<V>& c, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, const vector<K>& vdeg, K N, V P, V E, int L, int EF, K i, K n, vector<ThreadInfo*>& threads, FV fv) {
  if (EF!=LI_NORM) return 0;
  int MAX_THREADS = omp_get_max_threads();
  parallelBlockPthread(MAX_THREADS, [&](int t) {
    int& l = threads[t]->iteration;
    while (l<L) {
      V C0 = (1-P)/N;
      pagerankBarrierfreeCalculateRanksPthread(e, a, r, f, xv, xe, C0, E, i, n, threads[t], fv); ++l;  // update ranks of vertices
      if (!ASYNC) swap(a, r);                            // final ranks in (r)
      if (pagerankBarrierfreeConverged(e, i, n)) break;  // check tolerance
    }
    threads[t]->stop = timeNow();
  });
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
inline PagerankResult<V> pagerankBarrierfreePthread(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  using K = typename H::key_type;
  K    N  = xt.order();  if (N==0) return {};
  auto ks = vertexKeys(xt);
  return pagerankOmp<ASYNC>(xt, q, o, ks, 0, N, pagerankBarrierfreePthreadLoop<ASYNC, DEAD, K, V, FV>, fv);
}

/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @returns pagerank result
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V>
inline PagerankResult<V> pagerankBarrierfreePthread(const H& xt, const vector<V> *q, const PagerankOptions<V>& o) {
  auto fv = [](ThreadInfo *thread, auto v) {};
  return pagerankBarrierfreePthread<ASYNC, DEAD>(xt, q, o, fv);
}
