#ifndef PBRT_STD_H
#define PBRT_STD_H

#define __both__ __device__ __host__

namespace pbrt{

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
inline void SWAP(T& t1, T& t2) {
    T tmp(t1);
    t1=t2;
    t2=tmp;
}

__both__
inline int isinf(double x)
{
    union { uint64_t u; double f; } ieee754;
    ieee754.f = x;
    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) == 0x7ff00000 &&
           ( (unsigned)ieee754.u == 0 );
}

__both__
inline int isinf(float x)
{
    union { uint32_t u; float f; } ieee754;
    ieee754.f = x;
    return ( ieee754.u & 0x7fffffff == 0x7f800000 );
}

__both__
inline int isnan(double x)
{
    union { uint64_t u; double f; } ieee754;
    ieee754.f = x;
    return ( (unsigned)(ieee754.u >> 32) & 0x7fffffff ) +
           ( (unsigned)ieee754.u != 0 ) > 0x7ff00000;
}

__both__
inline int isnan(float x)
{
    union { uint32_t u; float f; } ieee754;
    ieee754.f = x;
    return ( ieee754.u & 0x7fffffff >  0x7f800000 );
}

template <class T> class numeric_limits{
  public:
    static constexpr T lowest();
    static constexpr T max();
};

template<> class numeric_limits<int>{
  public:
    __both__ static constexpr int lowest() { return -2147483648; }
    __both__ static constexpr int max() { return 2147483647; }
};

template<> class numeric_limits<float>{
  public:
    __both__ static constexpr float lowest() { return -3.40282e+38; }
    __both__ static constexpr float max() { return 3.40282e+38; }
};

template<> class numeric_limits<double>{
  public:
    __both__ static constexpr double lowest() { return -1.79769e+308; }
    __both__ static constexpr double max() { return 1.79769e+308; }
};

} // namespace pbrt
#endif
