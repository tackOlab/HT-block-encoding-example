#pragma once

#include <cstdint>
#include <cstdlib>

#define round_up(x, n) (((x) + (n)-1) & (-n))
#define round_down(x, n) ((x) & (-n))
#define ceil_int(a, b) ((a) + ((b)-1)) / (b)

#if defined(__arm64__) || defined(__arm__) || defined(__aarch64__)
  #include <arm_acle.h>
  #if defined(__ARM_NEON__)
    #include <arm_neon.h>
  #endif
#elif defined(_MSC_VER) || defined(__MINGW64__)
  #include <intrin.h>
#else
  #include <x86intrin.h>
#endif

static inline size_t popcount32(uintmax_t num) {
  size_t precision = 0;
#if defined(_MSC_VER)
  precision = __popcnt(static_cast<uint32_t>(num));
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
  precision = _popcnt32(num);
#else
  while (num != 0) {
    if (1 == (num & 1)) {
      precision++;
    }
    num >>= 1;
  }
#endif
  return precision;
}

static inline uint32_t int_log2(const uint32_t x) {
  uint32_t y;
#if defined(_MSC_VER)
  unsigned long tmp;
  _BitScanReverse(&tmp, x);
  y = tmp;
#else
  y         = 31 - __builtin_clz(x);
#endif
  return (x == 0) ? 0 : y;
}

static inline uint32_t count_leading_zeros(const uint32_t x) {
  uint32_t y;
#if defined(_MSC_VER)
  y = __lzcnt(x);
#elif defined(__AVX2__)
  y         = _lzcnt_u32(x);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  y      = __builtin_clz(x);
#elif defined(__ARM_FEATURE_CLZ)
  y = __builtin_clz(x);
#else
  y = 31 - int_log2(x);
#endif
  return (x == 0) ? 31 : y;
}

static inline void* aligned_mem_alloc(size_t size, size_t align) {
  void* result;
#if defined(__INTEL_COMPILER)
  result = _mm_malloc(size, align);
#elif defined(_MSC_VER)
  result    = _aligned_malloc(size, align);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  result = __mingw_aligned_malloc(size, align);
#else
  if (posix_memalign(&result, align, size)) {
    result = nullptr;
  }
#endif
  return result;
}

static inline void aligned_mem_free(void* ptr) {
#if defined(__INTEL_COMPILER)
  _mm_free(ptr);
#elif defined(_MSC_VER)
  _aligned_free(ptr);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  __mingw_aligned_free(ptr);
#else
  free(ptr);
#endif
}