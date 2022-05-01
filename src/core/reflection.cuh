
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

#ifndef PBRT_CORE_REFLECTION_H
#define PBRT_CORE_REFLECTION_H

// core/reflection.h*
#include "pbrt.cuh"
#include "geometry.cuh"
#include "shape.cuh"
#include "spectrum.cuh"

namespace pbrt {

// Reflection Declarations
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT);
Spectrum FrConductor(Float cosThetaI, const Spectrum &etaI,
                     const Spectrum &etaT, const Spectrum &k);

// BSDF Inline Functions
__both__
inline Float CosTheta(const Vector3f &w) { return w.z; }
__both__
inline Float Cos2Theta(const Vector3f &w) { return w.z * w.z; }
__both__
inline Float AbsCosTheta(const Vector3f &w) { return abs(w.z); }
__both__
inline Float Sin2Theta(const Vector3f &w) {
    return max((Float)0, (Float)1 - Cos2Theta(w));
}
__both__
inline Float SinTheta(const Vector3f &w) { return std::sqrt(Sin2Theta(w)); }
__both__
inline Float TanTheta(const Vector3f &w) { return SinTheta(w) / CosTheta(w); }
__both__
inline Float Tan2Theta(const Vector3f &w) {
    return Sin2Theta(w) / Cos2Theta(w);
}
__both__
inline Float CosPhi(const Vector3f &w) {
    Float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}
__both__
inline Float SinPhi(const Vector3f &w) {
    Float sinTheta = SinTheta(w);
    return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}
__both__
inline Float Cos2Phi(const Vector3f &w) { return CosPhi(w) * CosPhi(w); }
__both__
inline Float Sin2Phi(const Vector3f &w) { return SinPhi(w) * SinPhi(w); }
__both__
inline Float CosDPhi(const Vector3f &wa, const Vector3f &wb) {
    Float waxy = wa.x * wa.x + wa.y * wa.y;
    Float wbxy = wb.x * wb.x + wb.y * wb.y;
    if (waxy == 0 || wbxy == 0)
        return 1;
    return Clamp((wa.x * wb.x + wa.y * wb.y) / std::sqrt(waxy * wbxy), -1, 1);
}
__both__
inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) {
    return -wo + 2 * Dot(wo, n) * n;
}
__both__
inline bool Refract(const Vector3f &wi, const Normal3f &n, Float eta,
                    Vector3f *wt) {
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    Float cosThetaI = Dot(n, wi);
    Float sin2ThetaI = max(Float(0), Float(1 - cosThetaI * cosThetaI));
    Float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
    return true;
}
__both__
inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) {
    return w.z * wp.z > 0;
}
__both__
inline bool SameHemisphere(const Vector3f &w, const Normal3f &wp) {
    return w.z * wp.z > 0;
}

// BSDF Declarations
enum BxDFType {
    BSDF_REFLECTION = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE = 1 << 2,
    BSDF_GLOSSY = 1 << 3,
    BSDF_SPECULAR = 1 << 4,
    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
};

struct FourierBSDFTable {
    // FourierBSDFTable Public Data
    Float eta;
    int mMax;
    int nChannels;
    int nMu;
    Float *mu;
    int *m;
    int *aOffset;
    Float *a;
    Float *a0;
    Float *cdf;
    Float *recip;

    ~FourierBSDFTable() {
        delete[] mu;
        delete[] m;
        delete[] aOffset;
        delete[] a;
        delete[] a0;
        delete[] cdf;
        delete[] recip;
    }

    // FourierBSDFTable Public Methods
    static bool Read(const std::string &filename, FourierBSDFTable *table);
    const Float *GetAk(int offsetI, int offsetO, int *mptr) const {
        *mptr = m[offsetO * nMu + offsetI];
        return a + aOffset[offsetO * nMu + offsetI];
    }
    bool GetWeightsAndOffset(Float cosTheta, int *offset,
                             Float weights[4]) const;
};

class BSDF {
  public:
    // BSDF Public Methods
    __both__
    BSDF(const SurfaceInteraction &si, Float eta = 1)
        : eta(eta),
          ns(si.shading.n),
          ng(si.n),
          ss(Normalize(si.shading.dpdu)),
          ts(Cross(ns, ss)) {}
    __both__
    void Add(BxDF *b) {
        assert(nBxDFs < MaxBxDFs);
        bxdfs[nBxDFs++] = b;
    }
    __both__
    int NumComponents(BxDFType flags = BSDF_ALL) const;
    __both__
    Vector3f WorldToLocal(const Vector3f &v) const {
        return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    }
    __both__
    Vector3f LocalToWorld(const Vector3f &v) const {
        return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                        ss.y * v.x + ts.y * v.y + ns.y * v.z,
                        ss.z * v.x + ts.z * v.y + ns.z * v.z);
    }
    __both__
    Spectrum f(const Vector3f &woW, const Vector3f &wiW,
               BxDFType flags = BSDF_ALL) const;
    __both__
    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2,
                 BxDFType flags = BSDF_ALL) const;
    __both__
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples,
                 BxDFType flags = BSDF_ALL) const;
    __both__
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType type = BSDF_ALL,
                      BxDFType *sampledType = nullptr) const;
    __both__
    Float Pdf(const Vector3f &wo, const Vector3f &wi,
              BxDFType flags = BSDF_ALL) const;
    std::string ToString() const;

    // BSDF Public Data
    const Float eta;

  private:
    // BSDF Private Methods
    ~BSDF() {}

    // BSDF Private Data
    const Normal3f ns, ng;
    const Vector3f ss, ts;
    int nBxDFs = 0;
    static PBRT_CONSTEXPR int MaxBxDFs = 8;
    BxDF *bxdfs[MaxBxDFs];
    friend class MixMaterial;
};

inline std::ostream &operator<<(std::ostream &os, const BSDF &bsdf) {
    os << bsdf.ToString();
    return os;
}

// BxDF Declarations
class BxDF {
  public:
    // BxDF Interface
    ~BxDF() {}
    __both__
    BxDF(BxDFType type) : type(type) {}
    __both__
    bool MatchesFlags(BxDFType t) const { return (type & t) == type; }
    __both__
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi,
                              const Point2f &sample, Float *pdf,
                              BxDFType *sampledType = nullptr) const;
    __both__
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;

    // BxDF Public Data
    const BxDFType type;

    __both__
    BxDF(const Spectrum &R)
        : type(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
    __both__
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    __both__
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    __both__
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    // LambertianReflection Private Data
    const Spectrum R;
};

inline std::ostream &operator<<(std::ostream &os, const BxDF &bxdf) {
    os << bxdf.ToString();
    return os;
}

// class LambertianReflection : public BxDF {
//   public:
//     // LambertianReflection Public Methods
//     __both__
//     LambertianReflection(const Spectrum &R)
//         : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
//     __both__
//     Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
//     __both__
//     Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
//     __both__
//     Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
//     std::string ToString() const;

//   private:
//     // LambertianReflection Private Data
//     const Spectrum R;
// };

// BSDF Inline Method Definitions
__both__
inline int BSDF::NumComponents(BxDFType flags) const {
    int num = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) ++num;
    return num;
}

}  // namespace pbrt

#endif  // PBRT_CORE_REFLECTION_H
