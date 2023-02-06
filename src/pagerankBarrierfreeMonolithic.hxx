#pragma once
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "components.hxx"
#include "pagerank.hxx"
#include "pagerankBarrierfree.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;




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
inline PagerankResult<V> pagerankBarrierfreeMonolithicOmp(const G& x, const H& xt, const vector<V> *q, const PagerankOptions<V>& o, const PagerankData<G> *C, FV fv) {
  using K = typename H::key_type;
  K     N = xt.order();  if (N==0) return {};
  auto ks = joinValues(componentsD(x, xt, C));
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
inline PagerankResult<V> pagerankBarrierfreeMonolithicDynamicTraversalOmp(const G& x, const H& xt, const G& y, const H& yt, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions, const vector<V> *q, const PagerankOptions<V>& o, const PagerankData<G> *C, FV fv) {
  K    N    = yt.order();  if (N==0) return {};
  auto ks   = joinValues(componentsD(y, yt, C));
  auto vaff = compressContainer(y, pagerankAffectedTraversal(x, y, deletions, insertions), ks);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return pagerankOmp<ASYNC>(yt, q, o, ks, 0, N, pagerankBarrierfreeOmpLoop<ASYNC, DEAD, K, V, FV, decltype(fa)>, fv, fa);
}
