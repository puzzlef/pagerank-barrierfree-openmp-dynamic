#pragma once
#include <cstddef>
#include <utility>
#include <functional>
#include <algorithm>
#include <atomic>
#include <vector>
#include <pthread.h>

using std::function;
using std::atomic;
using std::pair;
using std::vector;
using std::make_pair;
using std::min;




// PARALLEL BLOCK
// --------------

inline void* parallelBlockHandlerPthread(void *data) {
  using DATA = pair<int, function<void(int)>>;
  const DATA& x = *((DATA*) data);
  x.second(x.first);
  return nullptr;
}


template <class F>
inline void parallelBlockPthread(int numThreads, F fn) {
  using  DATA = pair<int, function<void(int)>>;
  vector<DATA>      threadData(numThreads);
  vector<pthread_t> threads(numThreads);
  for (int i=0; i<numThreads; ++i) {
    threadData[i] = make_pair(i, function<void(int)>(fn));
    pthread_create(&threads[i], nullptr, parallelBlockHandlerPthread, &threadData[i]);
  }
  for (int i=0; i<numThreads; ++i)
    pthread_join(threads[i], nullptr);
}




// PARALLEL FOR
// ------------

template <class K, class F>
inline void parallelForDynamicNoWaitPthread(K begin, K end, F fn) {
  const  K CHUNK_SIZE = 2048;
  static atomic<K> xbegin(0), xend(0);
  K vbegin = 0, vend = 0;
  if (xbegin==0) xbegin.compare_exchange_weak(vbegin, begin);
  if (xend  ==0) xend  .compare_exchange_weak(vend,   end);
  for (K ibegin=xbegin; ibegin>=xend;) {
    K iend = min(ibegin + CHUNK_SIZE, K(xend));
    if (!xbegin.compare_exchange_weak(ibegin, iend)) continue;
    for (K i=ibegin; i<iend; ++i)
      fn(i);
  }
}
