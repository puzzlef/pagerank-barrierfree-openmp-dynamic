#pragma once
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"
#include "pagerankMonolithicSeq.hxx"

using std::vector;
using std::swap;
using std::min;




// PAGERANK-LOOP
// -------------

template <bool O, bool D, bool F, class T>
int pagerankBarrierfreeOmpLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<int>& vfrom, const vector<int>& efrom, const vector<int>& vdata, int i, int n, int N, T p, T E, int L, int EF) {
  float l = 0;
  if (!O) return 0;
  // Ordered approach
  int TS = omp_get_max_threads();
  int DN = ceilDiv(n, TS);
  #pragma omp parallel for schedule(static, 1) reduction(+:l)
  for (int t=0; t<TS; t++) {
    int  i1 = i+t*DN, I1 = min(i1+DN, i+n), n1 = I1-i1;
    int  l1 = pagerankMonolithicSeqLoopU<O, D, T, F>(a, r, c, f, vfrom, efrom, vdata, i1, n1, N, p, E, L, EF);
    l += l1 * n1/float(n);
  }
  return int(l + 0.5f);
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using multiple threads (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankBarrierfreeOmp(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = vertexKeys(xt);
  return pagerankOmp(xt, ks, 0, N, pagerankBarrierfreeOmpLoopU<O, D, F, T>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankBarrierfreeOmp(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  return pagerankBarrierfreeOmp<O, D, F>(x, xt, q, o);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankBarrierfreeOmpDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N = yt.order();                             if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = dynamicInVertices(x, xt, y, yt);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp(yt, ks, 0, n, pagerankBarrierfreeOmpLoopU<O, D, F, T>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankBarrierfreeOmpDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankBarrierfreeOmpDynamic<O, D, F>(x, xt, y, yt, q, o);
}
