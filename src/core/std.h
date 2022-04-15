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

} // namespace pbrt
#endif
