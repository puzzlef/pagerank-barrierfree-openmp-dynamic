#pragma once
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "transpose.hxx"
#include "components.hxx"
#include "sort.hxx"
#include "pagerank.hxx"
#include "pagerankBarrierfree.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;




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
 * @param ns vertex counts (per level)
 * @param threads information on each thread (updated)
 * @param fv per vertex processing (thread, vertex)
 * @param fa is vertex affected? (vertex)
 * @returns iterations performed
 */
template <bool ASYNC=false, bool DEAD=false, class K, class V, class NS, class FV, class FA>
inline int pagerankBarrierfreeLevelwiseOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, vector<V>& c, const vector<V>& f, const vector<size_t>& xv, const vector<K>& xe, const vector<K>& vdeg, K N, V P, V E, int L, int EF, K i, NS ns, vector<ThreadInfo*>& threads, FV fv, FA fa) {
  if (EF!=LI_NORM) return 0;
  pagerankBarrierfreeInitializeConvergedOmp(e, i, sumValues(ns), fa);
  #pragma omp parallel
  {
    K      j = i;
    int    t = omp_get_thread_num();
    float& l = threads[t]->iteration;
    for (auto n : ns) {
      while (l<L) {
        V C0 = DEAD? pagerankBarrierfreeTeleportOmp(r, vdeg, P, N) : (1-P)/N;
        pagerankBarrierfreeCalculateRanksOmp(e, a, r, f, xv, xe, C0, E, j, n, threads[t], fv, fa); l += float(n)/N;  // update ranks of vertices
        if (!ASYNC) swap(a, r);                            // final ranks in (r)
        if (pagerankBarrierfreeConverged(e, j, n)) break;  // check tolerance
        if (threads[t]->crashed) break;                    // simulate crash
      }
      j += n;
    }
    threads[t]->stop = timeNow();
  }
  int l = int(threadInfosMaxIteration(threads));
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
template <bool ASYNC=false, bool DEAD=false, class G, class H, class V, class FV>
inline PagerankResult<V> pagerankBarrierfreeLevelwiseOmp(const G& x, const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const PagerankData<G> *C, FV fv) {
  using K = typename H::key_type;
  K     N = xt.order();  if (N==0) return {};
  auto cs = componentsD(x, xt, C);
  auto  b = blockgraphD(x, cs, C);
  auto bt = blockgraphTransposeD(b, C);
  auto gs = levelwiseGroupedComponentsFrom(cs, b, bt);
  auto ks = joinValues(gs);
  auto fa = [](K u) { return true; };
  vector<K> ns;
  for (const auto& g : gs)
    ns.push_back(K(g.size()));
  return pagerankOmp<ASYNC>(xt, q, o, ks, 0, ns, pagerankBarrierfreeLevelwiseOmpLoop<ASYNC, DEAD, K, V, decltype(ns), FV, decltype(fa)>, fv, fa);
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
inline PagerankResult<V> pagerankBarrierfreeLevelwiseDynamicTraversalOmp(const G& x, const H& xt, const G& y, const H& yt, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions, const vector<V> *q, const PagerankOptions<V>& o, const PagerankData<G> *C, FV fv) {
  K    N    = yt.order();  if (N==0) return {};
  auto cs   = componentsD(y, yt, C);
  auto b    = blockgraphD(y, cs, C);
  auto bt   = blockgraphTransposeD(b, C);
  auto gs   = levelwiseGroupedComponentsFrom(cs, b, bt);
  auto ks   = joinValues(gs);
  auto vaff = compressContainer(y, pagerankAffectedTraversal(x, y, deletions, insertions), ks);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  vector<K> ns;
  for (const auto& g : gs)
    ns.push_back(K(g.size()));
  return pagerankOmp<ASYNC>(yt, q, o, ks, 0, ns, pagerankBarrierfreeLevelwiseOmpLoop<ASYNC, DEAD, K, V, decltype(ns), FV, decltype(fa)>, fv, fa);
}
