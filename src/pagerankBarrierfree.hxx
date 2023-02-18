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

#ifdef OPENMP
/**
 * Find total teleport contribution from each vertex (inc. dead ends).
 * @param xt transpose of original graph
 * @param r rank of each vertex
 * @param P damping factor [0.85]
 * @returns common teleport rank contribution to each vertex
 */
template <class H, class V>
inline V pagerankBarrierfreeTeleportOmp(const H& xt, const vector<V>& r, V P) {
  using  K = typename H::key_type;
  size_t S = xt.span();
  size_t N = xt.order();
  V a = (1-P)/N;
  #pragma omp for schedule(auto) reduction(+:a) nowait
  for (K u=0; u<S; ++u) {
    if (!xt.hasVertex(u)) continue;
    K   d = xt.vertexValue(u);
    if (d==0) a += P * r[u]/N;
  }
  return a;
}
#endif




// PAGERANK CALCULATE
// ------------------
// For rank calculation from in-edges.

/**
 * Calculate rank for a given vertex, and get the change in rank value.
 * @param a current rank of each vertex (output)
 * @param xt transpose of original graph
 * @param r previous rank of each vertex
 * @param f rank scaling factor for each vertex
 * @param v given vertex
 * @param C0 common teleport rank contribution to each vertex
 * @param P damping factor [0.85]
 * @returns change between previous and current rank value
 */
template <class H, class K, class V>
inline V pagerankCalculateRankDelta(vector<V>& a, const H& xt, const vector<V>& r, K v, V C0, V P) {
  V av = C0, rv = r[v];
  xt.forEachEdgeKey(v, [&](auto u) {
    K d = xt.vertexData(u);
    av += P * r[u]/d;
  });
  a[v] = av;
  return abs(av - rv);
}


#ifdef OPENMP
/**
 * Calculate ranks for vertices in a graph.
 * @param e change in rank for each vertex below tolerance? (updated)
 * @param a current rank of each vertex (updated)
 * @param xt transpose of original graph
 * @param r previous rank of each vertex
 * @param C0 common teleport rank contribution to each vertex
 * @param P damping factor [0.85]
 * @param E tolerance [10^-10]
 * @param thread information on current thread (updated)
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 */
template <class H, class V, class FV, class FA>
inline void pagerankBarrierfreeCalculateRanksOmp(vector<int>& e, vector<V>& a, const H& xt, const vector<V>& r, V C0, V P, V E, ThreadInfo *thread, FV fv, FA fa) {
  using  K = typename H::key_type;
  size_t S = xt.span();
  #pragma omp for schedule(dynamic, 2048) nowait
  for (K v=0; v<S; ++v) {
    if (!xt.hasVertex(v) || !fa(v)) continue;
    V ev = pagerankCalculateRankDelta(a, xt, r, v, C0, P);
    if (ev<=E && e[v]==0) e[v] = 1;  // LI_NORM
    fv(thread, v);
  }
}
#endif




// PAGERANK CONVERGED
// ------------------
// For convergence check.

#ifdef OPENMP
/**
 * Mark unaffected vertices as converged.
 * @param e change in rank for each vertex below tolerance?
 * @param xt transpose of original graph
 * @param fa is vertex affected? (vertex)
 */
template <class H, class FA>
inline void pagerankBarrierfreeInitializeConvergedOmp(vector<int>& e, const H& xt, FA fa) {
  using  K = typename H::key_type;
  size_t S = xt.span();
  #pragma omp for schedule(auto) nowait
  for (K u=0; u<S; ++u)
    if (!xt.hasVertex(v) || !fa(u)) e[u] = 1;
}
#endif


/**
 * Check if ranks of all vertices have converged.
 * @param e change in rank for each vertex below tolerance?
 * @param xt transpose of original graph
 * @returns ranks converged?
 */
template <class H>
inline bool pagerankBarrierfreeConverged(const vector<int>& e, const H& xt) {
  using  K = typename H::key_type;
  size_t S = xt.span();
  for (K u=0; u<S; ++u)
    if (xt.hasVertex(u) && !e[u]) return false;
  return true;
}




// PAGERANK AFFECTED (TRAVERSAL)
// -----------------------------

#ifdef OPENMP
/**
 * Find affected vertices due to a batch update.
 * @param vis affected flags (output)
 * @param x original graph
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 */
template <class B, class G, class K>
inline void pagerankBarrierfreeAffectedTraversalOmp(vector<B>& vis, const G& x, const G& y, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions) {
  auto fn = [](K u) {};
  #pragma omp for schedule(auto) nowait
  for (size_t i=0, I=deletions.size(); i<I; ++i) {
    K u = get<0>(deletions[i]);
    dfsVisitedForEachW(vis, x, u, fn);
  }
  #pragma omp for schedule(auto) nowait
  for (size_t i=0, I=insertions.size(); i<I; ++i) {
    K u = get<0>(insertions[i]);
    dfsVisitedForEachW(vis, y, u, fn);
  }
}
#endif




// PAGERANK LOOP
// -------------

#ifdef OPENMP
/**
 * Perform PageRank iterations upon a graph.
 * @param e change in rank for each vertex below tolerance? (unused)
 * @param a current rank of each vertex (updated)
 * @param r previous rank of each vertex (updated)
 * @param xt transpose of original graph
 * @param P damping factor [0.85]
 * @param E tolerance [10^-10]
 * @param L max. iterations [500]
 * @param EF error function (L1/L2/LI)
 * @param threads information on each thread (updated)
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 * @param fp preprocessing to perform
 * @returns iterations performed
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V, class FV, class FA, class FP>
inline int pagerankBarrierfreeOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads, FV fv, FA fa, FP fp) {
  if (EF!=LI_NORM) return 0;
  #pragma omp parallel
  {
    fp();
    pagerankBarrierfreeInitializeConvergedOmp(e, xt, fa);
    int  t = omp_get_thread_num();
    int& l = threads[t]->iteration;
    while (l<L) {
      V C0 = DEAD? pagerankBarrierfreeTeleportOmp(xt, r, P) : (1-P)/N;
      pagerankBarrierfreeCalculateRanksOmp(e, a, xt, r, C0, P, E, threads[t], fv, fa); ++l;  // update ranks of vertices
      if (!ASYNC) swap(a, r);                          // final ranks in (r)
      if (pagerankBarrierfreeConverged(e, xt)) break;  // check tolerance
      if (threads[t]->crashed) break;                  // simulate crash
    }
    threads[t]->stop = timeNow();
  }
  int l = threadInfosMaxIteration(threads);
  if (!ASYNC && (l & 1)==1) swap(a, r);
  return l;
}
#endif




// STATIC/NAIVE-DYNAMIC PAGERANK
// -----------------------------

#ifdef OPENMP
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
  if  (xt.empty()) return {};
  return pagerankOmp<ASYNC>(xt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [](K u) { return true; };
    auto fp = []() {};
    return pagerankBarrierfreeOmpLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa, fp);
  });
}
#endif




// TRAVERSAL-BASED DYNAMIC PAGERANK
// --------------------------------

#ifdef OPENMP
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
  if (xt.empty()) return {};
  vector<char> vaff(max(x.span(), y.span()));
  return pagerankOmp<ASYNC>(yt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [&](K u) { return vaff[u]==1; };
    auto fp = [&]()    { pagerankBarrierfreeAffectedTraversalOmpW(vaff, x, y, deletions, insertions); };
    return pagerankBasicOmpLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa, fp);
  });
}
#endif
