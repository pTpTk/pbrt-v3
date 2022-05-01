
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

// core/reflection.cpp*
#include "reflection.cuh"
#include "spectrum.cuh"
#include "sampler.cuh"
#include "sampling.cuh"
#include "interpolation.cuh"
#include "scene.cuh"
#include "interaction.cuh"
#include "stats.cuh"
#include <stdarg.h>

namespace pbrt {

// BxDF Utility Functions
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    // Potentially swap indices of refraction
    bool entering = cosThetaI > 0.f;
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    Float sinThetaI = pbrt::math::sqrt(std::max((Float)0, 1 - cosThetaI * cosThetaI));
    Float sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    Float cosThetaT = pbrt::math::sqrt(std::max((Float)0, 1 - sinThetaT * sinThetaT));
    Float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    Float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
Spectrum FrConductor(Float cosThetaI, const Spectrum &etai,
                     const Spectrum &etat, const Spectrum &k) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    Spectrum eta = etat / etai;
    Spectrum etak = k / etai;

    Float cosThetaI2 = cosThetaI * cosThetaI;
    Float sinThetaI2 = 1. - cosThetaI2;
    Spectrum eta2 = eta * eta;
    Spectrum etak2 = etak * etak;

    Spectrum t0 = eta2 - etak2 - sinThetaI2;
    Spectrum a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
    Spectrum t1 = a2plusb2 + cosThetaI2;
    Spectrum a = Sqrt(0.5f * (a2plusb2 + t0));
    Spectrum t2 = (Float)2 * cosThetaI * a;
    Spectrum Rs = (t1 - t2) / (t1 + t2);

    Spectrum t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    Spectrum t4 = t2 * sinThetaI2;
    Spectrum Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
}

// // BxDF Method Definitions
// Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const {
//     return R * InvPi;
// }

std::string BxDF::ToString() const {
    return std::string("[ LambertianReflection R: ") + R.ToString() +
           std::string(" ]");
}

// Spectrum BxDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
//                         Float *pdf, BxDFType *sampledType) const {
//     // Cosine-sample the hemisphere, flipping the direction if necessary
//     *wi = CosineSampleHemisphere(u);
//     if (wo.z < 0) wi->z *= -1;
//     *pdf = Pdf(wo, *wi);
//     return f(wo, *wi);
// }

// Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
//     return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
// }

// Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
//     Spectrum r(0.);
//     for (int i = 0; i < nSamples; ++i) {
//         // Estimate one term of $\rho_\roman{hd}$
//         Vector3f wi;
//         Float pdf = 0;
//         Spectrum f = Sample_f(w, &wi, u[i], &pdf);
//         if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
//     }
//     return r / nSamples;
// }

// Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
//     Spectrum r(0.f);
//     for (int i = 0; i < nSamples; ++i) {
//         // Estimate one term of $\rho_\roman{hh}$
//         Vector3f wo, wi;
//         wo = UniformSampleHemisphere(u1[i]);
//         Float pdfo = UniformHemispherePdf(), pdfi = 0;
//         Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
//         if (pdfi > 0)
//             r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
//     }
//     return r / (Pi * nSamples);
// }

// BSDF Method Definitions
// Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW,
//                  BxDFType flags) const {
//     // ProfilePhase pp(Prof::BSDFEvaluation);
//     Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
//     if (wo.z == 0) return 0.;
//     bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
//     Spectrum f(0.f);
//     for (int i = 0; i < nBxDFs; ++i)
//         if (bxdfs[i]->MatchesFlags(flags) &&
//             ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
//              (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
//             f += bxdfs[i]->f(wo, wi);
//     return f;
// }

Spectrum BSDF::rho(int nSamples, const Point2f *samples1,
                   const Point2f *samples2, BxDFType flags) const {
    Spectrum ret(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(nSamples, samples1, samples2);
    return ret;
}

Spectrum BSDF::rho(const Vector3f &woWorld, int nSamples, const Point2f *samples,
                   BxDFType flags) const {
    Vector3f wo = WorldToLocal(woWorld);
    Spectrum ret(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(wo, nSamples, samples);
    return ret;
}

// Spectrum BSDF::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
//                         const Point2f &u, Float *pdf, BxDFType type,
//                         BxDFType *sampledType) const {
//     // ProfilePhase pp(Prof::BSDFSampling);
//     // Choose which _BxDF_ to sample
//     int matchingComps = NumComponents(type);
//     if (matchingComps == 0) {
//         *pdf = 0;
//         if (sampledType) *sampledType = BxDFType(0);
//         return Spectrum(0);
//     }
//     int comp =
//         min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

//     // Get _BxDF_ pointer for chosen component
//     BxDF *bxdf = nullptr;
//     int count = comp;
//     for (int i = 0; i < nBxDFs; ++i)
//         if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
//             bxdf = bxdfs[i];
//             break;
//         }
//     // CHECK(bxdf != nullptr);
//     // VLOG(2) << "BSDF::Sample_f chose comp = " << comp << " / matching = " <<
//         // matchingComps << ", bxdf: " << bxdf->ToString();

//     // Remap _BxDF_ sample _u_ to $[0,1)^2$
//     Point2f uRemapped(min(u[0] * matchingComps - comp, OneMinusEpsilon),
//                       u[1]);

//     // Sample chosen _BxDF_
//     Vector3f wi, wo = WorldToLocal(woWorld);
//     if (wo.z == 0) return 0.;
//     *pdf = 0;
//     if (sampledType) *sampledType = bxdf->type;
//     Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
//     // VLOG(2) << "For wo = " << wo << ", sampled f = " << f << ", pdf = "
//     //         << *pdf << ", ratio = " << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.))
//     //         << ", wi = " << wi;
//     if (*pdf == 0) {
//         if (sampledType) *sampledType = BxDFType(0);
//         return 0;
//     }
//     *wiWorld = LocalToWorld(wi);

//     // Compute overall PDF with all matching _BxDF_s
//     if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
//         for (int i = 0; i < nBxDFs; ++i)
//             if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type))
//                 *pdf += bxdfs[i]->Pdf(wo, wi);
//     if (matchingComps > 1) *pdf /= matchingComps;

//     // Compute value of BSDF for sampled direction
//     if (!(bxdf->type & BSDF_SPECULAR)) {
//         bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
//         f = 0.;
//         for (int i = 0; i < nBxDFs; ++i)
//             if (bxdfs[i]->MatchesFlags(type) &&
//                 ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
//                  (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
//                 f += bxdfs[i]->f(wo, wi);
//     }
//     // VLOG(2) << "Overall f = " << f << ", pdf = " << *pdf << ", ratio = "
//     //         << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.));
//     return f;
// }

// Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
//                 BxDFType flags) const {
//     // ProfilePhase pp(Prof::BSDFPdf);
//     if (nBxDFs == 0.f) return 0.f;
//     Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
//     if (wo.z == 0) return 0.;
//     Float pdf = 0.f;
//     int matchingComps = 0;
//     for (int i = 0; i < nBxDFs; ++i)
//         if (bxdfs[i]->MatchesFlags(flags)) {
//             ++matchingComps;
//             pdf += bxdfs[i]->Pdf(wo, wi);
//         }
//     Float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
//     return v;
// }

std::string BSDF::ToString() const {
    std::string s = StringPrintf("[ BSDF eta: %f nBxDFs: %d", eta, nBxDFs);
    for (int i = 0; i < nBxDFs; ++i)
        s += StringPrintf("\n  bxdfs[%d]: ", i) + bxdfs[i]->ToString();
    return s + std::string(" ]");
}

}  // namespace pbrt
