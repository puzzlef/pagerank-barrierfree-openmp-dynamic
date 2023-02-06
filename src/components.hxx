#pragma once
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "dfs.hxx"

using std::vector;




// COMPONENTS
// ----------
// Finds Strongly Connected Components (SCC) using Kosaraju's algorithm.

template <class G, class H>
inline auto components(const G& x, const H& xt) {
  using K = typename G::key_type;
  vector2d<K> a; vector<K> vs;
  // original dfs
  auto vis = createContainer(x, bool());
  x.forEachVertexKey([&](auto u) {
    if (vis[u]) return;
    dfsEndVisitedForEachW(vis, x, u, [&](auto v) { vs.push_back(v); });
  });
  // transpose dfs
  fillValueU(vis, false);
  while (!vs.empty()) {
    auto u = vs.back(); vs.pop_back();
    if (vis[u]) continue;
    vector<K> c;
    dfsVisitedForEachW(vis, xt, u, [&](auto v) { c.push_back(v); });
    a.push_back(move(c));
  }
  return a;
}
