#pragma once
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "pagerank.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::swap;




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
inline V pagerankBarrierfreeTeleportOmp(const vector<V>& r, const vector<K>& vdeg, V P, K N) {
  V a = (1-P)/N;
  #pragma omp for schedule(auto) reduction(+:a) nowait
  for (K u=0; u<N; ++u)
    if (vdeg[u]==0) a += P * r[u]/N;
  return a;
}




// PAGERANK CALCULATE
// ------------------
// For rank calculation from in-edges.

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
 * @param fa is vertex affected? (vertex)
 */
template <class K, class V, class FV, class FA>
inline void pagerankBarrierfreeCalculateRanksOmp(vector<int>& e, vector<V>& a, const vector<V>& r, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, V C0, V E, K i, K n, ThreadInfo *thread, FV fv, FA fa) {
  #pragma omp for schedule(dynamic, 2048) nowait
  for (K v=i; v<i+n; ++v) {
    if (!fa(v)) continue;
    V ev = pagerankCalculateRankDelta(a, r, f, xv, xe, v, C0);
    if (ev<=E && e[v]==0) e[v] = 1;  // LI_NORM
    fv(thread, v);
  }
}




// PAGERANK CONVERGED
// ------------------
// For convergence check.

/**
 * Mark unaffected vertices as converged.
 * @param e change in rank for each vertex below tolerance?
 * @param i vertex start
 * @param n vertex count
 * @param fa is vertex affected? (vertex)
 */
template <class K, class FA>
inline void pagerankBarrierfreeInitializeConvergedOmp(vector<int>& e, K i, K n, FA fa) {
  #pragma omp parallel for schedule(auto)
  for (K u=i; u<i+n; ++u)
    if (!fa(u)) e[u] = 1;
}


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




// PAGERANK AFFECTED (TRAVERSAL)
// -----------------------------

#ifdef OPENMP
/**
 * Find affected vertices due to a batch update.
 * @param x original graph
 * @param y updated graph
 * @param ft is vertex affected? (u)
 * @returns affected flags
 */
template <bool B=char, class G, class K>
inline auto pagerankAffectedTraversalBarrierfreeOmp(const G& x, const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions) {
  auto fn = [](K u) {};
  vector<B> vis(max(x.span(), y.span()));
  #pragma omp for schedule(auto)
  for (size_t i=0, I=deletions.size(); i<I; ++i) {
    K u = get<0>(deletions[i]);
    dfsVisitedForEachW(vis, x, u, fn);
  }
  #pragma omp for schedule(auto)
  for (size_t i=0, I=insertions.size(); i<I; ++i) {
    K u = get<0>(insertions[i]);
    dfsVisitedForEachW(vis, y, u, fn);
  }
  return vis;
}
#endif




// PAGERANK LOOP
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
 * @param fa is vertex affected? (vertex)
 * @returns iterations performed
 */
template <bool ASYNC=false, bool DEAD=false, class K, class V, class FV, class FA>
inline int pagerankBarrierfreeOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, vector<V>& c, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, const vector<K>& vdeg, K N, V P, V E, int L, int EF, K i, K n, vector<ThreadInfo*>& threads, FV fv, FA fa) {
  if (EF!=LI_NORM) return 0;
  pagerankBarrierfreeInitializeConvergedOmp(e, i, n, fa);
  #pragma omp parallel
  {
    int  t = omp_get_thread_num();
    int& l = threads[t]->iteration;
    while (l<L) {
      V C0 = DEAD? pagerankBarrierfreeTeleportOmp(r, vdeg, P, N) : (1-P)/N;
      pagerankBarrierfreeCalculateRanksOmp(e, a, r, f, xv, xe, C0, E, i, n, threads[t], fv, fa); ++l;  // update ranks of vertices
      if (!ASYNC) swap(a, r);                            // final ranks in (r)
      if (pagerankBarrierfreeConverged(e, i, n)) break;  // check tolerance
      if (threads[t]->crashed) break;                    // simulate crash
    }
    threads[t]->stop = timeNow();
  }
  int l = threadInfosMaxIteration(threads);
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
  K     N = xt.order();  if (N==0) return {};
  auto ks = vertexKeys(xt);
  auto fa = [](K u) { return true; };
  return pagerankOmp<ASYNC>(xt, q, o, ks, 0, N, pagerankBarrierfreeOmpLoop<ASYNC, DEAD, K, V, FV, decltype(fa)>, fv, fa);
}




// TRAVERSAL-BASED DYNAMIC PAGERANK
// --------------------------------

/**
 * Find the rank of each vertex in a dynamic graph.
 * @param x original graph
 * @param xt transpose of original graph
 * @param y updated graph
 * @param yt transpose of updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial ranks
 * @param o pagerank options
 * @param fv per vertex processing (thread, vertex)
 * @returns pagerank result
 */
template <bool ASYNC=false, bool DEAD=false, class G, class H, class K, class V, class FV>
inline PagerankResult<V> pagerankBarrierfreeDynamicTraversalOmp(const G& x, const H& xt, const G& y, const H& yt, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  K    N    = yt.order();  if (N==0) return {};
  auto ks   = vertexKeys(yt);
  auto vaff = compressContainer(y, pagerankAffectedVerticesTraversalBarrierfreeOmp(x, deletions, insertions), ks);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return pagerankOmp<ASYNC>(yt, q, o, ks, 0, N, pagerankBarrierfreeOmpLoop<ASYNC, DEAD, K, V, FV, decltype(fa)>, fv, fa);
}
