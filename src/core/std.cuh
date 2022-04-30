#ifndef PBRT_STD_H
#define PBRT_STD_H

namespace pbrt{
#include <utility>

#define __both__ __device__ __host__

namespace math{
__both__
inline float sqrt(float x) {
   #ifdef __CUDA_ARCH__
      return sqrtf(x);
   #else
      return std::sqrt(x);
   #endif
}

__both__
inline double sqrt(double x) {
   #ifdef __CUDA_ARCH__
      return sqrt(x);
   #else
      return std::sqrt(x);
   #endif
}

__both__
inline long double sqrt(long double x) {
   #ifdef __CUDA_ARCH__
      return sqrt(x);
   #else
      return std::sqrt(x);
   #endif
}

__both__
inline long double abs(long double x) {
   #ifdef __CUDA_ARCH__
      return fabsf(x);
   #else
      return std::abs(x);
   #endif
}

__both__
inline float abs(float x) {
   #ifdef __CUDA_ARCH__
      return fabsf(x);
   #else
      return std::abs(x);
   #endif
}

__both__
inline double abs(double x) {
   #ifdef __CUDA_ARCH__
      return fabs(x);
   #else
      return std::abs(x);
   #endif
}

__both__
inline int isnan(float x){
   #ifdef __CUDA_ARCH__
      return isnan(x);
   #else
      return std::isnan(x);
   #endif
}

__both__
inline int isnan(double x){
   #ifdef __CUDA_ARCH__
      return isnan(x);
   #else
      return std::isnan(x);
   #endif
}

__both__
inline bool isnan(int x){return false;}

}  // namespace math

template<class T>
__both__
inline const T& min(const T& a, const T& b)
{
   return (b < a) ? b : a;
}

template<class T>
__both__
inline const T& max(const T& a, const T& b)
{
   return (b > a) ? b : a;
}

template<typename T> 
__both__
inline void Swap(T& t1, T& t2) {
   T tmp(t1);
   t1=t2;
   t2=tmp;
}

__both__
inline bool isinf(double x)
{
   union { uint64_t u; double f; } ieee754;
   ieee754.f = x;
   return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) == 0x7ff00000 &&
          ( (unsigned)ieee754.u == 0 );
}

__both__
inline bool isinf(float x)
{
   union { uint32_t u; float f; } ieee754;
   ieee754.f = x;
   return ( ieee754.u & 0x7fffffff == 0x7f800000 );
}

__both__
inline bool isinf(int x){return false;}

// __both__
// inline bool isnan(double x)
// {
//    union { uint64_t u; double f; } ieee754;
//    ieee754.f = x;
//    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) +
//           ( (unsigned)ieee754.u != 0 ) > 0x7ff00000;
// }

// __both__
// inline bool isnan(float x)
// {
//    union { uint32_t u; float f; } ieee754;
//    ieee754.f = x;
//    return ( ieee754.u & 0x7fffffff >  0x7f800000 );
// }



template <class T> class numeric_limits{
  public:
    __both__ static constexpr T lowest();
    __both__ static constexpr T max();
    __both__ static constexpr T epsilon();
    __both__ static constexpr T infinity();
};

template<> class numeric_limits<int>{
  public:
    __both__ static constexpr int lowest() { return -2147483648; }
    __both__ static constexpr int max() { return 2147483647; }
    __both__ static constexpr int epsilon() { return 0; }
    __both__ static constexpr int infinity() { return 0; }
};

template<> class numeric_limits<float>{
  public:
    __both__ static constexpr float lowest() { return -3.40282e+38; }
    __both__ static constexpr float max() { return 3.40282e+38; }
    __both__ static constexpr float epsilon() { return 1.19209e-07; }
    __both__ static constexpr float infinity() { return max() + 2; }
};

template<> class numeric_limits<double>{
  public:
    __both__ static constexpr double lowest() { return -1.79769e+308; }
    __both__ static constexpr double max() { return 1.79769e+308; }
    __both__ static constexpr double epsilon() { return 2.22045e-16; }
    __both__ static constexpr double infinity() { return max() + 2; }
};


template<typename T>
class shared_ptr {
public:
   __both__
   shared_ptr() : ptr ((T*)malloc(sizeof(T))) {}
   __both__
   shared_ptr(T *_ptr) : ptr ((T*)malloc(sizeof(T))) {
      *ptr = *_ptr;
      //ptr = new (ptr) T(*_ptr); // placement new doesn't work for virtual classes
      // delete _ptr; ?? Maybe delete the sink pointer
   }

   // std::shared_ptr to pbrt::gpu::shared_ptr needs to be invoked from host
   __host__
   shared_ptr(const std::shared_ptr<T>& _ptr) : ptr (_ptr.get()) {}
   __both__
   T& operator*() const noexcept { return *(get()); }
   __both__
   T* operator->() const noexcept { return get(); }
   __both__
   T *get() const noexcept { return ptr; }
   __both__
   operator bool() const noexcept { return get() != nullptr; }
private:
   T *ptr;
};

template<typename T>
class unique_ptr {
public:
   __both__
   unique_ptr() : ptr ((T*)malloc(sizeof(T))) {}
   __both__
   unique_ptr(T *_ptr) : ptr ((T*)malloc(sizeof(T))) {
      *ptr = *_ptr;
      //ptr = new (ptr) T(*_ptr); // placement new doesn't work for virtual classes
      // delete _ptr; ?? Maybe delete the sink pointer
   }

   // std::shared_ptr to pbrt::gpu::shared_ptr needs to be invoked from host
   __host__
   unique_ptr(const std::unique_ptr<T>& _ptr) : ptr (_ptr.get()) {}
   __both__
   T& operator*() const noexcept { return *(get()); }
   __both__
   T* operator->() const noexcept { return get(); }
   __both__
   T *get() const noexcept { return ptr; }
   __both__
   operator bool() const noexcept { return get() != nullptr; }
private:
   T *ptr;
};

template<typename T>
__both__
shared_ptr<T> make_shared(T obj) {
   shared_ptr<T> ptr;
   *ptr = obj;
   return ptr;
}

namespace utils{

template<typename T>
__both__
int get_buffer_size(T const * const ptr) {
  int count = 0;
  for (auto it = ptr; it != nullptr; ++it)
    ++count;
  return count;
}

template<typename T>
__both__
void swap(T& x, T& y) {
  T temp = std::move(x);
  x = std::move(y);
  y = std::move(temp);
}

}  // namespace utils

}  // namespace pbrt
#endif