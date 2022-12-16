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
inline void parallelForDynamicNoWaitPthread(K begin, K end, K chunkSize, F fn) {
  static atomic<K>      xbegin(0), xend(0);
  static atomic<size_t> xround(0);
  thread_local  size_t  tround(0);
  if  (++tround<xround) return;
  for (K obegin=xbegin, oend=xend; obegin>=oend && tround>xround;) {
    xend    .compare_exchange_strong(oend,   min(begin, obegin));
    xbegin  .compare_exchange_strong(obegin, begin);
    if (xend.compare_exchange_strong(oend,   end)) ++xround;
    break;
  }
  for (K ibegin=xbegin; ibegin<xend;) {
    K iend = min(ibegin + chunkSize, K(xend));
    if (!xbegin.compare_exchange_strong(ibegin, iend)) continue;
    for (K i=ibegin; i<iend; ++i)
      fn(i);
  }
}
