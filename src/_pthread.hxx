#pragma once
#include <cstddef>
#include <algorithm>
#include <atomic>
#include <vector>
#include <pthread.h>

using std::atomic;
using std::vector;
using std::min;




// PARALLEL BLOCK
// --------------

template <class F>
inline void parallelBlockHandler(void *data) {
  int threadNum = *((int*) data);
  F fn; fn(threadNum);
}


template <class F>
inline void parallelBlock(int numThreads, F fn) {
  vector<int>       threadNums(numThreads);
  vector<pthread_t> threads   (numThreads);
  for (int i=0; i<numThreads; ++i) {
    threadNums[i] = i;
    pthread_create(&threads[i], nullptr, parallelBlockHandler<F>, &threadNums[i]);
  }
  for (int i=0; i<numThreads; ++i)
    pthread_join(threads[i]);
}




// PARALLEL FOR
// ------------

template <class K, class F>
inline void parallelForDynamicNoWaitPthread(K begin, K end, F fn) {
  const  K CHUNK_SIZE = 2048;
  static atomic<K> xbegin(0), xend(0);
  if (xbegin==0) xbegin.compare_exchange_weak(0, begin);
  if (xend  ==0) xend  .compare_exchange_weak(0, end);
  while (true) {
    K ibegin = xbegin;
    if (ibegin >= xend) break;
    K iend   = min(ibegin + CHUNK_SIZE, xend);
    if (!xbegin.compare_exchange_weak(ibegin, iend)) continue;
    for (K i=ibegin; i<iend; ++i)
      fn(i);
  }
}
