#pragma once
#include <utility>
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "pagerank.hxx"

using std::tuple;
using std::vector;
using std::swap;
using std::max;




// PAGERANK LOOP
// -------------

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
 * @returns iterations performed
 */
template <bool ASYNC=false, bool DEAD=false, class H, class V, class FV, class FA>
inline int pagerankBasicSeqLoop(vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads, FV fv, FA fa) {
  using  K = typename H::key_type;
  size_t N = xt.order();
  int l = 0;
  while (l<L) {
    V C0 = DEAD? pagerankTeleport(xt, r, P) : (1-P)/N;
    pagerankCalculateRanks(a, xt, r, C0, P, threads[0], fv, fa); ++l;  // update ranks of vertices
    V el = pagerankError(a, r, EF);  // compare previous and current ranks
    if (!ASYNC) swap(a, r);          // final ranks in (r)
    if (el<E) break;                 // check tolerance
    if (threads[0]->crashed) break;  // simulate crash
  }
  return l;
}


#ifdef OPENMP
template <bool ASYNC=false, bool DEAD=false, class H, class V, class FV, class FA>
inline int pagerankBasicOmpLoop(vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads, FV fv, FA fa) {
  using  K = typename H::key_type;
  size_t N = xt.order();
  int l = 0;
  while (l<L) {
    V C0 = DEAD? pagerankTeleportOmp(xt, r, P) : (1-P)/N;
    pagerankCalculateRanksOmp(a, xt, r, C0, P, threads, fv, fa); ++l;  // update ranks of vertices
    V el = pagerankErrorOmp(a, r, EF);  // compare previous and current ranks
    if (!ASYNC) swap(a, r);             // final ranks in (r)
    if (el<E) break;                    // check tolerance
    if (threadInfosCrashedCount(threads) > 0) break;  // simulate crash
  }
  return l;
}
#endif




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
inline PagerankResult<V> pagerankBasicSeq(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  using K = typename H::key_type;
  if  (xt.empty()) return {};
  return pagerankSeq<ASYNC>(xt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [](K u) { return true; };
    return pagerankBasicSeqLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa);
  });
}


#ifdef OPENMP
template <bool ASYNC=false, bool DEAD=false, class H, class V, class FV>
inline PagerankResult<V> pagerankBasicOmp(const H& xt, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  using K = typename H::key_type;
  if  (xt.empty()) return {};
  return pagerankOmp<ASYNC>(xt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [](K u) { return true; };
    return pagerankBasicOmpLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa);
  });
}
#endif




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
inline PagerankResult<V> pagerankBasicDynamicTraversalSeq(const G& x, const H& xt, const G& y, const H& yt, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  if (xt.empty()) return {};
  vector<bool> vaff(max(x.span(), y.span()));
  return pagerankSeq<ASYNC>(yt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [&](K u) { return vaff[u]==true; };
    pagerankAffectedTraversalW(vaff, x, y, deletions, insertions);
    return pagerankBasicSeqLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa);
  });
}


#ifdef OPENMP
template <bool ASYNC=false, bool DEAD=false, class G, class H, class K, class V, class FV>
inline PagerankResult<V> pagerankBasicDynamicTraversalOmp(const G& x, const H& xt, const G& y, const H& yt, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K>>& insertions, const vector<V> *q, const PagerankOptions<V>& o, FV fv) {
  if (xt.empty()) return {};
  vector<char> vaff(max(x.span(), y.span()));
  return pagerankOmp<ASYNC>(yt, q, o, [&](vector<int>& e, vector<V>& a, vector<V>& r, const H& xt, V P, V E, int L, int EF, vector<ThreadInfo*>& threads) {
    auto fa = [&](K u) { return vaff[u]==1; };
    pagerankAffectedTraversalOmpW(vaff, x, y, deletions, insertions);
    return pagerankBasicOmpLoop<ASYNC, DEAD>(e, a, r, xt, P, E, L, EF, threads, fv, fa);
  });
}
#endif
