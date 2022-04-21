
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_GEOMETRY_H
#define PBRT_CORE_GEOMETRY_H

// core/geometry.h*
#include "pbrt.cuh"
#include "stringprint.cuh"
#include <iterator>

namespace pbrt {
namespace gpu {

template <typename T>
__both__
inline bool isNaN(const T x) {
    return isnan(x);
}
template <>
__both__
inline bool isNaN(const int x) {
    return false;
}

// Vector Declarations
template <typename T>
class Vector2 {
  public:
    // Vector2 Public Methods
    __both__
    Vector2() { x = y = 0; }
    __both__
    Vector2(T xx, T yy) : x(xx), y(yy) { assert(!HasNaNs()); }
    __both__
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }
    __both__
    explicit Vector2(const Point2<T> &p);
    __both__
    explicit Vector2(const Point3<T> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    __both__
    Vector2(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
    }
    __both__
    Vector2<T> &operator=(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        return *this;
    }
#endif  // !NDEBUG
    __both__
    Vector2<T> operator+(const Vector2<T> &v) const {
        assert(!v.HasNaNs());
        return Vector2(x + v.x, y + v.y);
    }
    __both__
    Vector2<T> &operator+=(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    __both__
    Vector2<T> operator-(const Vector2<T> &v) const {
        assert(!v.HasNaNs());
        return Vector2(x - v.x, y - v.y);
    }
    __both__
    Vector2<T> &operator-=(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    __both__
    bool operator==(const Vector2<T> &v) const { return x == v.x && y == v.y; }
    __both__
    bool operator!=(const Vector2<T> &v) const { return x != v.x || y != v.y; }
    template <typename U>
    __both__
    Vector2<T> operator*(U f) const {
        return Vector2<T>(f * x, f * y);
    }

    template <typename U>
    __both__
    Vector2<T> &operator*=(U f) {
        assert(!isNaN(f));
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    __both__
    Vector2<T> operator/(U f) const {
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Vector2<T>(x * inv, y * inv);
    }

    template <typename U>
    __both__
    Vector2<T> &operator/=(U f) {
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    __both__
    Vector2<T> operator-() const { return Vector2<T>(-x, -y); }
    __both__
    T operator[](int i) const {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    __both__
    T &operator[](int i) {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    __both__
    Float LengthSquared() const { return x * x + y * y; }
    __both__
    Float Length() const { return std::sqrt(LengthSquared()); }

    // Vector2 Public Data
    T x, y;
};

template <typename T>
class Vector3 {
  public:
    // Vector3 Public Methods
    __both__
    T operator[](int i) const {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    __both__
    T &operator[](int i) {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    __both__
    Vector3() { x = y = z = 0; }
    __both__
    Vector3(T x, T y, T z) : x(x), y(y), z(z) { assert(!HasNaNs()); }
    __both__
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    __both__
    explicit Vector3(const Point3<T> &p);
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    __both__
    Vector3(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
    }
    __both__
    Vector3<T> &operator=(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }
#endif  // !NDEBUG
    __both__
    Vector3<T> operator+(const Vector3<T> &v) const {
        assert(!v.HasNaNs());
        return Vector3(x + v.x, y + v.y, z + v.z);
    }
    __both__
    Vector3<T> &operator+=(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    __both__
    Vector3<T> operator-(const Vector3<T> &v) const {
        assert(!v.HasNaNs());
        return Vector3(x - v.x, y - v.y, z - v.z);
    }
    __both__
    Vector3<T> &operator-=(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    __both__
    bool operator==(const Vector3<T> &v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    __both__
    bool operator!=(const Vector3<T> &v) const {
        return x != v.x || y != v.y || z != v.z;
    }
    template <typename U>
    __both__
    Vector3<T> operator*(U s) const {
        return Vector3<T>(s * x, s * y, s * z);
    }
    template <typename U>
    __both__
    Vector3<T> &operator*=(U s) {
        assert(!isNaN(s));
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    template <typename U>
    __both__
    Vector3<T> operator/(U f) const {
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Vector3<T>(x * inv, y * inv, z * inv);
    }
    template <typename U>
    __both__
    Vector3<T> &operator/=(U f) {
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    __both__
    Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }
    __both__
    Float LengthSquared() const { return x * x + y * y + z * z; }
    __both__
    Float Length() const { return std::sqrt(LengthSquared()); }
    __both__
    explicit Vector3(const Normal3<T> &n);

    // Vector3 Public Data
    T x, y, z;
};

typedef Vector2<Float> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector3<Float> Vector3f;
typedef Vector3<int> Vector3i;

// Point Declarations
template <typename T>
class Point2 {
  public:
    // Point2 Public Methods
    __both__
    explicit Point2(const Point3<T> &p) : x(p.x), y(p.y) { assert(!HasNaNs()); }
    __both__
    Point2() { x = y = 0; }
    __both__
    Point2(T xx, T yy) : x(xx), y(yy) { assert(!HasNaNs()); }

    template <typename U>
    __both__
    explicit Point2(const Point2<U> &p) {
        x = (T)p.x;
        y = (T)p.y;
        assert(!HasNaNs());
    }

    template <typename U>
    __both__
    explicit Point2(const Vector2<U> &p) {
        x = (T)p.x;
        y = (T)p.y;
        assert(!HasNaNs());
    }

    template <typename U>
    __both__
    explicit operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }

#ifndef NDEBUG
    __both__
    Point2(const Point2<T> &p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
    }
    __both__
    Point2<T> &operator=(const Point2<T> &p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        return *this;
    }
#endif  // !NDEBUG
    __both__
    Point2<T> operator+(const Vector2<T> &v) const {
        assert(!v.HasNaNs());
        return Point2<T>(x + v.x, y + v.y);
    }
    __both__
    Point2<T> &operator+=(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    __both__
    Vector2<T> operator-(const Point2<T> &p) const {
        assert(!p.HasNaNs());
        return Vector2<T>(x - p.x, y - p.y);
    }
    __both__
    Point2<T> operator-(const Vector2<T> &v) const {
        assert(!v.HasNaNs());
        return Point2<T>(x - v.x, y - v.y);
    }
    __both__
    Point2<T> operator-() const { return Point2<T>(-x, -y); }
    __both__
    Point2<T> &operator-=(const Vector2<T> &v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    __both__
    Point2<T> &operator+=(const Point2<T> &p) {
        assert(!p.HasNaNs());
        x += p.x;
        y += p.y;
        return *this;
    }
    __both__
    Point2<T> operator+(const Point2<T> &p) const {
        assert(!p.HasNaNs());
        return Point2<T>(x + p.x, y + p.y);
    }
    template <typename U>
    __both__
    Point2<T> operator*(U f) const {
        return Point2<T>(f * x, f * y);
    }
    template <typename U>
    __both__
    Point2<T> &operator*=(U f) {
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    __both__
    Point2<T> operator/(U f) const {
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Point2<T>(inv * x, inv * y);
    }
    template <typename U>
    __both__
    Point2<T> &operator/=(U f) {
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    __both__
    T operator[](int i) const {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    __both__
    T &operator[](int i) {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    __both__
    bool operator==(const Point2<T> &p) const { return x == p.x && y == p.y; }
    __both__
    bool operator!=(const Point2<T> &p) const { return x != p.x || y != p.y; }
    __both__
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }

    // Point2 Public Data
    T x, y;
};

template <typename T>
class Point3 {
  public:
    // Point3 Public Methods
    __both__
    Point3() { x = y = z = 0; }
    __both__
    Point3(T x, T y, T z) : x(x), y(y), z(z) { assert(!HasNaNs()); }
    template <typename U>
    __both__
    explicit Point3(const Point3<U> &p)
        : x((T)p.x), y((T)p.y), z((T)p.z) {
        assert(!HasNaNs());
    }
    template <typename U>
    __both__
    explicit operator Vector3<U>() const {
        return Vector3<U>(x, y, z);
    }
#ifndef NDEBUG
    __both__
    Point3(const Point3<T> &p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }
    __both__
    Point3<T> &operator=(const Point3<T> &p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
#endif  // !NDEBUG
    __both__
    Point3<T> operator+(const Vector3<T> &v) const {
        assert(!v.HasNaNs());
        return Point3<T>(x + v.x, y + v.y, z + v.z);
    }
    __both__
    Point3<T> &operator+=(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    __both__
    Vector3<T> operator-(const Point3<T> &p) const {
        assert(!p.HasNaNs());
        return Vector3<T>(x - p.x, y - p.y, z - p.z);
    }
    __both__
    Point3<T> operator-(const Vector3<T> &v) const {
        assert(!v.HasNaNs());
        return Point3<T>(x - v.x, y - v.y, z - v.z);
    }
    __both__
    Point3<T> &operator-=(const Vector3<T> &v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    __both__
    Point3<T> &operator+=(const Point3<T> &p) {
        assert(!p.HasNaNs());
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    __both__
    Point3<T> operator+(const Point3<T> &p) const {
        assert(!p.HasNaNs());
        return Point3<T>(x + p.x, y + p.y, z + p.z);
    }
    template <typename U>
    __both__
    Point3<T> operator*(U f) const {
        return Point3<T>(f * x, f * y, f * z);
    }
    template <typename U>
    __both__
    Point3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    __both__
    Point3<T> operator/(U f) const {
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Point3<T>(inv * x, inv * y, inv * z);
    }
    template <typename U>
    __both__
    Point3<T> &operator/=(U f) {
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    __both__
    T operator[](int i) const {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    __both__
    T &operator[](int i) {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    __both__
    bool operator==(const Point3<T> &p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    __both__
    bool operator!=(const Point3<T> &p) const {
        return x != p.x || y != p.y || z != p.z;
    }
    __both__
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    __both__
    Point3<T> operator-() const { return Point3<T>(-x, -y, -z); }

    // Point3 Public Data
    T x, y, z;
};

typedef Point2<Float> Point2f;
typedef Point2<int> Point2i;
typedef Point3<Float> Point3f;
typedef Point3<int> Point3i;

// Normal Declarations
template <typename T>
class Normal3 {
  public:
    // Normal3 Public Methods
    __both__
    Normal3() { x = y = z = 0; }
    __both__
    Normal3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) { assert(!HasNaNs()); }
    __both__
    Normal3<T> operator-() const { return Normal3(-x, -y, -z); }
    __both__
    Normal3<T> operator+(const Normal3<T> &n) const {
        assert(!n.HasNaNs());
        return Normal3<T>(x + n.x, y + n.y, z + n.z);
    }
    __both__
    Normal3<T> &operator+=(const Normal3<T> &n) {
        assert(!n.HasNaNs());
        x += n.x;
        y += n.y;
        z += n.z;
        return *this;
    }
    __both__
    Normal3<T> operator-(const Normal3<T> &n) const {
        assert(!n.HasNaNs());
        return Normal3<T>(x - n.x, y - n.y, z - n.z);
    }
    __both__
    Normal3<T> &operator-=(const Normal3<T> &n) {
        assert(!n.HasNaNs());
        x -= n.x;
        y -= n.y;
        z -= n.z;
        return *this;
    }
    __both__
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    template <typename U>
    __both__
    Normal3<T> operator*(U f) const {
        return Normal3<T>(f * x, f * y, f * z);
    }

    template <typename U>
    __both__
    Normal3<T> &operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    __both__
    Normal3<T> operator/(U f) const {
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Normal3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    __both__
    Normal3<T> &operator/=(U f) {
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    __both__
    Float LengthSquared() const { return x * x + y * y + z * z; }
    __both__
    Float Length() const { return std::sqrt(LengthSquared()); }

#ifndef NDEBUG
    __both__
    Normal3<T>(const Normal3<T> &n) {
        assert(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
    }
    __both__
    Normal3<T> &operator=(const Normal3<T> &n) {
        assert(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
        return *this;
    }
#endif  // !NDEBUG
    __both__
    explicit Normal3<T>(const Vector3<T> &v) : x(v.x), y(v.y), z(v.z) {
        assert(!v.HasNaNs());
    }
    __both__
    bool operator==(const Normal3<T> &n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    __both__
    bool operator!=(const Normal3<T> &n) const {
        return x != n.x || y != n.y || z != n.z;
    }
    __both__
    T operator[](int i) const {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    __both__
    T &operator[](int i) {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    // Normal3 Public Data
    T x, y, z;
};

typedef Normal3<Float> Normal3f;

// Bounds Declarations
template <typename T>
class Bounds2 {
  public:
    // Bounds2 Public Methods
    __both__
    Bounds2() {
        T minNum = numeric_limits<T>::lowest();
        T maxNum = numeric_limits<T>::max();
        pMin = Point2<T>(maxNum, maxNum);
        pMax = Point2<T>(minNum, minNum);
    }
    __both__
    explicit Bounds2(const Point2<T> &p) : pMin(p), pMax(p) {}
    __both__
    Bounds2(const Point2<T> &p1, const Point2<T> &p2) {
        pMin = Point2<T>(min(p1.x, p2.x), min(p1.y, p2.y));
        pMax = Point2<T>(max(p1.x, p2.x), max(p1.y, p2.y));
    }
    template <typename U>
    __both__
    explicit operator Bounds2<U>() const {
        return Bounds2<U>((Point2<U>)pMin, (Point2<U>)pMax);
    }
    __both__
    Vector2<T> Diagonal() const { return pMax - pMin; }
    __both__
    T Area() const {
        Vector2<T> d = pMax - pMin;
        return (d.x * d.y);
    }
    __both__
    int MaximumExtent() const {
        Vector2<T> diag = Diagonal();
        if (diag.x > diag.y)
            return 0;
        else
            return 1;
    }
    __both__
    inline const Point2<T> &operator[](int i) const {
        assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }
    __both__
    inline Point2<T> &operator[](int i) {
        assert(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }
    __both__
    bool operator==(const Bounds2<T> &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    __both__
    bool operator!=(const Bounds2<T> &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    __both__
    Point2<T> Lerp(const Point2f &t) const {
        return Point2<T>(pbrt::gpu::Lerp(t.x, pMin.x, pMax.x),
                         pbrt::gpu::Lerp(t.y, pMin.y, pMax.y));
    }
    __both__
    Vector2<T> Offset(const Point2<T> &p) const {
        Vector2<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        return o;
    }
    __both__
    void BoundingSphere(Point2<T> *c, Float *rad) const {
        *c = (pMin + pMax) / 2;
        *rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
    }

    // Bounds2 Public Data
    Point2<T> pMin, pMax;
};

template <typename T>
class Bounds3 {
  public:
    // Bounds3 Public Methods
    __both__
    Bounds3() {
        T minNum = numeric_limits<T>::lowest();
        T maxNum = numeric_limits<T>::max();
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }
    __both__
    explicit Bounds3(const Point3<T> &p) : pMin(p), pMax(p) {}
    __both__
    Bounds3(const Point3<T> &p1, const Point3<T> &p2)
        : pMin(min(p1.x, p2.x), min(p1.y, p2.y),
               min(p1.z, p2.z)),
          pMax(max(p1.x, p2.x), max(p1.y, p2.y),
               max(p1.z, p2.z)) {}
    __both__
    const Point3<T> &operator[](int i) const;
    __both__
    Point3<T> &operator[](int i);
    __both__
    bool operator==(const Bounds3<T> &b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    __both__
    bool operator!=(const Bounds3<T> &b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    __both__
    Point3<T> Corner(int corner) const {
        assert(corner >= 0 && corner < 8);
        return Point3<T>((*this)[(corner & 1)].x,
                         (*this)[(corner & 2) ? 1 : 0].y,
                         (*this)[(corner & 4) ? 1 : 0].z);
    }
    __both__
    Vector3<T> Diagonal() const { return pMax - pMin; }
    __both__
    T SurfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    __both__
    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }
    __both__
    int MaximumExtent() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }
    __both__
    Point3<T> Lerp(const Point3f &t) const {
        return Point3<T>(pbrt::gpu::Lerp(t.x, pMin.x, pMax.x),
                         pbrt::gpu::Lerp(t.y, pMin.y, pMax.y),
                         pbrt::gpu::Lerp(t.z, pMin.z, pMax.z));
    }
    __both__
    Vector3<T> Offset(const Point3<T> &p) const {
        Vector3<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
        return o;
    }
    __both__
    void BoundingSphere(Point3<T> *center, Float *radius) const {
        *center = (pMin + pMax) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }
    template <typename U>
    __both__
    explicit operator Bounds3<U>() const {
        return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
    }
    __both__
    bool IntersectP(const gpu::Ray &ray, Float *hitt0 = nullptr,
                    Float *hitt1 = nullptr) const;
    __both__
    inline bool IntersectP(const gpu::Ray &ray, const Vector3f &invDir,
                           const int dirIsNeg[3]) const;

    // Bounds3 Public Data
    Point3<T> pMin, pMax;
};

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

class Bounds2iIterator : public std::forward_iterator_tag {
  public:
    __both__
    Bounds2iIterator(const Bounds2i &b, const Point2i &pt)
        : p(pt), bounds(&b) {}
    __both__
    Bounds2iIterator operator++() {
        advance();
        return *this;
    }
    __both__
    Bounds2iIterator operator++(int) {
        Bounds2iIterator old = *this;
        advance();
        return old;
    }
    __both__
    bool operator==(const Bounds2iIterator &bi) const {
        return p == bi.p && bounds == bi.bounds;
    }
    __both__
    bool operator!=(const Bounds2iIterator &bi) const {
        return p != bi.p || bounds != bi.bounds;
    }
    __both__
    Point2i operator*() const { return p; }

  private:
    __both__
    void advance() {
        ++p.x;
        if (p.x == bounds->pMax.x) {
            p.x = bounds->pMin.x;
            ++p.y;
        }
    }
    Point2i p;
    const Bounds2i *bounds;
};

// Ray Declarations
class Ray {
  public:
    // Ray Public Methods
    __both__
    Ray() : tMax(Infinity), time(0.f), medium(nullptr) {}
    __both__
    Ray(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
        Float time = 0.f, const Medium *medium = nullptr)
        : o(o), d(d), tMax(tMax), time(time), medium(medium) {}
    __both__
    Point3f operator()(Float t) const { return o + d * t; }
    __both__
    bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }

    // Ray Public Data
    Point3f o;
    Vector3f d;
    mutable Float tMax;
    Float time;
    const Medium *medium;
};

class RayDifferential : public Ray {
  public:
    // RayDifferential Public Methods
    __both__
    RayDifferential() { hasDifferentials = false; }
    __both__
    RayDifferential(const Point3f &o, const Vector3f &d, Float tMax = Infinity,
                    Float time = 0.f, const Medium *medium = nullptr)
        : Ray(o, d, tMax, time, medium) {
        hasDifferentials = false;
    }
    __both__
    RayDifferential(const Ray &ray) : Ray(ray) { hasDifferentials = false; }
    __both__
    bool HasNaNs() const {
        return Ray::HasNaNs() ||
               (hasDifferentials &&
                (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                 rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }
    __both__
    void ScaleDifferentials(Float s) {
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    Point3f rxOrigin, ryOrigin;
    Vector3f rxDirection, ryDirection;
};

// Geometry Inline Functions
template <typename T>
__both__
inline Vector3<T>::Vector3(const Point3<T> &p)
    : x(p.x), y(p.y), z(p.z) {
    assert(!HasNaNs());
}

template <typename T, typename U>
__both__
inline Vector3<T> operator*(U s, const Vector3<T> &v) {
    return v * s;
}
template <typename T>
__both__
Vector3<T> Abs(const Vector3<T> &v) {
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
__both__
inline T Dot(const Vector3<T> &v1, const Vector3<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
__both__
inline T AbsDot(const Vector3<T> &v1, const Vector3<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
__both__
inline Vector3<T> Cross(const Vector3<T> &v1, const Vector3<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
__both__
inline Vector3<T> Cross(const Vector3<T> &v1, const Normal3<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
__both__
inline Vector3<T> Cross(const Normal3<T> &v1, const Vector3<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
                      (v1x * v2y) - (v1y * v2x));
}

template <typename T>
__both__
inline Vector3<T> Normalize(const Vector3<T> &v) {
    return v / v.Length();
}
template <typename T>
__both__
T MinComponent(const Vector3<T> &v) {
    return min(v.x, min(v.y, v.z));
}

template <typename T>
__both__
T MaxComponent(const Vector3<T> &v) {
    return max(v.x, max(v.y, v.z));
}

template <typename T>
__both__
int MaxDimension(const Vector3<T> &v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template <typename T>
__both__
Vector3<T> Min(const Vector3<T> &p1, const Vector3<T> &p2) {
    return Vector3<T>(min(p1.x, p2.x), min(p1.y, p2.y),
                      min(p1.z, p2.z));
}

template <typename T>
__both__
Vector3<T> Max(const Vector3<T> &p1, const Vector3<T> &p2) {
    return Vector3<T>(max(p1.x, p2.x), max(p1.y, p2.y),
                      max(p1.z, p2.z));
}

template <typename T>
__both__
Vector3<T> Permute(const Vector3<T> &v, int x, int y, int z) {
    return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
__both__
inline void CoordinateSystem(const Vector3<T> &v1, Vector3<T> *v2,
                             Vector3<T> *v3) {
    if (std::abs(v1.x) > std::abs(v1.y))
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
    *v3 = Cross(v1, *v2);
}

template <typename T>
__both__
Vector2<T>::Vector2(const Point2<T> &p)
    : x(p.x), y(p.y) {
    assert(!HasNaNs());
}

template <typename T>
__both__
Vector2<T>::Vector2(const Point3<T> &p)
    : x(p.x), y(p.y) {
    assert(!HasNaNs());
}

template <typename T, typename U>
__both__
inline Vector2<T> operator*(U f, const Vector2<T> &v) {
    return v * f;
}
template <typename T>
__both__
inline Float Dot(const Vector2<T> &v1, const Vector2<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
__both__
inline Float AbsDot(const Vector2<T> &v1, const Vector2<T> &v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
__both__
inline Vector2<T> Normalize(const Vector2<T> &v) {
    return v / v.Length();
}
template <typename T>
__both__
Vector2<T> Abs(const Vector2<T> &v) {
    return Vector2<T>(std::abs(v.x), std::abs(v.y));
}

template <typename T>
__both__
inline Float Distance(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1 - p2).Length();
}

template <typename T>
__both__
inline Float DistanceSquared(const Point3<T> &p1, const Point3<T> &p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
__both__
inline Point3<T> operator*(U f, const Point3<T> &p) {
    assert(!p.HasNaNs());
    return p * f;
}

template <typename T>
__both__
Point3<T> Lerp(Float t, const Point3<T> &p0, const Point3<T> &p1) {
    return (1 - t) * p0 + t * p1;
}

template <typename T>
__both__
Point3<T> Min(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(min(p1.x, p2.x), min(p1.y, p2.y),
                     min(p1.z, p2.z));
}

template <typename T>
__both__
Point3<T> Max(const Point3<T> &p1, const Point3<T> &p2) {
    return Point3<T>(max(p1.x, p2.x), max(p1.y, p2.y),
                     max(p1.z, p2.z));
}

template <typename T>
__both__
Point3<T> Floor(const Point3<T> &p) {
    return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
__both__
Point3<T> Ceil(const Point3<T> &p) {
    return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
__both__
Point3<T> Abs(const Point3<T> &p) {
    return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
__both__
inline Float Distance(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1 - p2).Length();
}

template <typename T>
__both__
inline Float DistanceSquared(const Point2<T> &p1, const Point2<T> &p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
__both__
inline Point2<T> operator*(U f, const Point2<T> &p) {
    assert(!p.HasNaNs());
    return p * f;
}

template <typename T>
__both__
Point2<T> Floor(const Point2<T> &p) {
    return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
__both__
Point2<T> Ceil(const Point2<T> &p) {
    return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
__both__
Point2<T> Lerp(Float t, const Point2<T> &v0, const Point2<T> &v1) {
    return (1 - t) * v0 + t * v1;
}

template <typename T>
__both__
Point2<T> Min(const Point2<T> &pa, const Point2<T> &pb) {
    return Point2<T>(min(pa.x, pb.x), min(pa.y, pb.y));
}

template <typename T>
__both__
Point2<T> Max(const Point2<T> &pa, const Point2<T> &pb) {
    return Point2<T>(max(pa.x, pb.x), max(pa.y, pb.y));
}

template <typename T>
__both__
Point3<T> Permute(const Point3<T> &p, int x, int y, int z) {
    return Point3<T>(p[x], p[y], p[z]);
}

template <typename T, typename U>
__both__
inline Normal3<T> operator*(U f, const Normal3<T> &n) {
    return Normal3<T>(f * n.x, f * n.y, f * n.z);
}

template <typename T>
__both__
inline Normal3<T> Normalize(const Normal3<T> &n) {
    return n / n.Length();
}

template <typename T>
__both__
inline Vector3<T>::Vector3(const Normal3<T> &n)
    : x(n.x), y(n.y), z(n.z) {
    assert(!n.HasNaNs());
}

template <typename T>
__both__
inline T Dot(const Normal3<T> &n1, const Vector3<T> &v2) {
    assert(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

template <typename T>
__both__
inline T Dot(const Vector3<T> &v1, const Normal3<T> &n2) {
    assert(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

template <typename T>
__both__
inline T Dot(const Normal3<T> &n1, const Normal3<T> &n2) {
    assert(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T>
__both__
inline T AbsDot(const Normal3<T> &n1, const Vector3<T> &v2) {
    assert(!n1.HasNaNs() && !v2.HasNaNs());
    return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}

template <typename T>
__both__
inline T AbsDot(const Vector3<T> &v1, const Normal3<T> &n2) {
    assert(!v1.HasNaNs() && !n2.HasNaNs());
    return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}

template <typename T>
__both__
inline T AbsDot(const Normal3<T> &n1, const Normal3<T> &n2) {
    assert(!n1.HasNaNs() && !n2.HasNaNs());
    return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

template <typename T>
__both__
inline Normal3<T> Faceforward(const Normal3<T> &n, const Vector3<T> &v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T>
__both__
inline Normal3<T> Faceforward(const Normal3<T> &n, const Normal3<T> &n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
__both__
inline Vector3<T> Faceforward(const Vector3<T> &v, const Vector3<T> &v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
__both__
inline Vector3<T> Faceforward(const Vector3<T> &v, const Normal3<T> &n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template <typename T>
__both__
Normal3<T> Abs(const Normal3<T> &v) {
    return Normal3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
__both__
inline const Point3<T> &Bounds3<T>::operator[](int i) const {
    assert(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}

template <typename T>
__both__
inline Point3<T> &Bounds3<T>::operator[](int i) {
    assert(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
}

template <typename T>
__both__
Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
    Bounds3<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
__both__
Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    Bounds3<T> ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
__both__
Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    Bounds3<T> ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
__both__
bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}

template <typename T>
__both__
bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
            p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}

template <typename T>
__both__
bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
            p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z);
}

template <typename T, typename U>
__both__
inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
                      b.pMax + Vector3<T>(delta, delta, delta));
}

// Minimum squared distance from point to box; returns zero if point is
// inside.
template <typename T, typename U>
__both__
inline Float DistanceSquared(const Point3<T> &p, const Bounds3<U> &b) {
    Float dx = max({Float(0), b.pMin.x - p.x, p.x - b.pMax.x});
    Float dy = max({Float(0), b.pMin.y - p.y, p.y - b.pMax.y});
    Float dz = max({Float(0), b.pMin.z - p.z, p.z - b.pMax.z});
    return dx * dx + dy * dy + dz * dz;
}

template <typename T, typename U>
__both__
inline Float Distance(const Point3<T> &p, const Bounds3<U> &b) {
    return std::sqrt(DistanceSquared(p, b));
}
__both__
inline Bounds2iIterator begin(const Bounds2i &b) {
    return Bounds2iIterator(b, b.pMin);
}
__both__
inline Bounds2iIterator end(const Bounds2i &b) {
    // Normally, the ending point is at the minimum x value and one past
    // the last valid y value.
    Point2i pEnd(b.pMin.x, b.pMax.y);
    // However, if the bounds are degenerate, override the end point to
    // equal the start point so that any attempt to iterate over the bounds
    // exits out immediately.
    if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
        pEnd = b.pMin;
    return Bounds2iIterator(b, pEnd);
}

template <typename T>
__both__
Bounds2<T> Union(const Bounds2<T> &b, const Point2<T> &p) {
    Bounds2<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
__both__
Bounds2<T> Union(const Bounds2<T> &b, const Bounds2<T> &b2) {
    Bounds2<T> ret;
    ret.pMin = Min(b.pMin, b2.pMin);
    ret.pMax = Max(b.pMax, b2.pMax);
    return ret;
}

template <typename T>
__both__
Bounds2<T> Intersect(const Bounds2<T> &b1, const Bounds2<T> &b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    Bounds2<T> ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
__both__
bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb) {
    bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
    bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
    return (x && y);
}

template <typename T>
__both__
bool Inside(const Point2<T> &pt, const Bounds2<T> &b) {
    return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
            pt.y <= b.pMax.y);
}

template <typename T>
__both__
bool InsideExclusive(const Point2<T> &pt, const Bounds2<T> &b) {
    return (pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
            pt.y < b.pMax.y);
}

template <typename T, typename U>
__both__
Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
    return Bounds2<T>(b.pMin - Vector2<T>(delta, delta),
                      b.pMax + Vector2<T>(delta, delta));
}

template <typename T>
__both__
inline bool Bounds3<T>::IntersectP(const gpu::Ray &ray, Float *hitt0,
                                   Float *hitt1) const {
    Float t0 = 0, t1 = ray.tMax;
    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
        Float invRayDir = 1 / ray.d[i];
        Float tNear = (pMin[i] - ray.o[i]) * invRayDir;
        Float tFar = (pMax[i] - ray.o[i]) * invRayDir;

        // Update parametric interval from slab intersection $t$ values
        if (tNear > tFar) pbrt::gpu::Swap(tNear, tFar);

        // Update _tFar_ to ensure robust ray--bounds intersection
        tFar *= 1 + 2 * gamma(3);
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1) return false;
    }
    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}

template <typename T>
__both__
inline bool Bounds3<T>::IntersectP(const gpu::Ray &ray, const Vector3f &invDir,
                                   const int dirIsNeg[3]) const {
    const Bounds3f &bounds = *this;
    // Check for ray intersection against $x$ and $y$ slabs
    Float tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
    Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

    // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
    tMax *= 1 + 2 * gamma(3);
    tyMax *= 1 + 2 * gamma(3);
    if (tMin > tyMax || tyMin > tMax) return false;
    if (tyMin > tMin) tMin = tyMin;
    if (tyMax < tMax) tMax = tyMax;

    // Check for ray intersection against $z$ slab
    Float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
    Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

    // Update _tzMax_ to ensure robust bounds intersection
    tzMax *= 1 + 2 * gamma(3);
    if (tMin > tzMax || tzMin > tMax) return false;
    if (tzMin > tMin) tMin = tzMin;
    if (tzMax < tMax) tMax = tzMax;
    return (tMin < ray.tMax) && (tMax > 0);
}
__both__
inline Point3f OffsetRayOrigin(const Point3f &p, const Vector3f &pError,
                               const Normal3f &n, const Vector3f &w) {
    Float d = Dot(Abs(n), pError);
    Vector3f offset = d * Vector3f(n);
    if (Dot(w, n) < 0) offset = -offset;
    Point3f po = p + offset;
    // Round offset point _po_ away from _p_
    for (int i = 0; i < 3; ++i) {
        if (offset[i] > 0)
            po[i] = NextFloatUp(po[i]);
        else if (offset[i] < 0)
            po[i] = NextFloatDown(po[i]);
    }
    return po;
}
__both__
inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                    cosTheta);
}
__both__
inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
                                   const Vector3f &x, const Vector3f &y,
                                   const Vector3f &z) {
    return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
           cosTheta * z;
}
__both__
inline Float SphericalTheta(const Vector3f &v) {
    return std::acos(Clamp(v.z, -1, 1));
}
__both__
inline Float SphericalPhi(const Vector3f &v) {
    Float p = std::atan2(v.y, v.x);
    return (p < 0) ? (p + 2 * Pi) : p;
}

}   // namespace gpu
}  // namespace pbrt

#endif  // PBRT_CORE_GEOMETRY_H
