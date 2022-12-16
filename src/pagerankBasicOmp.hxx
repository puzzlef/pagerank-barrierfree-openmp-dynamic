#pragma once
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "transpose.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"

using std::vector;
using std::swap;




// PAGERANK-LOOP
// -------------

/**
 * Perform PageRank iterations upon a graph.
 * @param e change in rank for each vertex below tolerance? (unused)
 * @param a current rank of each vertex (updated)
 * @param r previous rank of each vertex (updated)
 * @param c rank contribution from each vertex (updated)
 * @param f rank scaling factor for each vertex
 * @param xv edge offsets for each vertex in the graph
 * @param xe target vertices for each edge in the graph
 * @param vdeg out-degree of each vertex
 * @param N total number of vertices
 * @param P damping factor [0.85]
 * @param E tolerance [10^-10]
 * @param L max. iterations [500]
 * @param EF error function (L1/L2/LI)
 * @param i vertex start
 * @param n vertex count
 * @param threads information on each thread (updated)
 * @param fv per vertex processing (thread, vertex)
 * @returns iterations performed
 */
template <bool ASYNC=false, bool DEAD=false, class K, class V, class FV>
int pagerankBasicOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, vector<V>& c, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, const vector<K>& vdeg, K N, V P, V E, int L, int EF, K i, K n, vector<ThreadInfo*>& threads, FV fv) {
  int l = 0;
  while (l<L) {
    V C0 = DEAD? pagerankTeleportOmp(r, vdeg, P, N) : (1-P)/N;
    pagerankCalculateRanksOmp(a, c, xv, xe, C0, i, n, threads, fv); ++l;  // update ranks of vertices
    multiplyValuesOmpW(c, a, f, i, n);        // update partial contributions (c)
    V el = pagerankErrorOmp(a, r, EF, i, n);  // compare previous and current ranks
    if (!ASYNC) swap(a, r);                   // final ranks in (r)
    if (el<E) break;                          // check tolerance
  }
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
inline PagerankResult<V> pagerankBasicOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  using K = typename H::key_type;
  K    N  = xt.order();  if (N==0) return {};
  auto ks = vertexKeys(xt);
  return pagerankOmp<ASYNC>(xt, q, o, ks, 0, N, pagerankBasicOmpLoop<ASYNC, DEAD, K, V, FV>, fv);
}

/**
 * Find the rank of each vertex in a graph.
 * @param xt transpose of original graph
 * @param q initial ranks
 * @param o pagerank options
 * @returns pagerank result
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V>
inline PagerankResult<V> pagerankBasicOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o) {
  auto fv = [](ThreadInfo *thread, auto v) {};
  return pagerankBasicOmp<ASYNC, DEAD>(xt, q, o, fv);
}
