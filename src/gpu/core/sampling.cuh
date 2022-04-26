
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

#ifndef PBRT_CORE_SAMPLING_H
#define PBRT_CORE_SAMPLING_H

// core/sampling.h*
#include "pbrt.cuh"
#include "geometry.cuh"
#include "rng.cuh"
#include <algorithm>

namespace pbrt {
namespace gpu {
__both__
inline int FindIntervalSampling(int size, const Float *nodes, const Float x) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (nodes[middle] <= x) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return Clamp(first - 1, 0, size - 2);
}

// Sampling Declarations
__both__
void StratifiedSample1D(Float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
__both__
void StratifiedSample2D(Point2f *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
void LatinHypercube(Float *samples, int nSamples, int nDim, RNG &rng);
struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const Float *f, int n) {
        // Compute integral of step function at $x_i$

        // only need to realloc memory if the desired size is different from last call
        // spoiler: it never is, a vector was used in place of a statically sized array for no reason
        if(N != n+1){
            N = n + 1;
            if(cdf != nullptr) free(cdf);
            cdf = (Float*) malloc(N * sizeof(Float));

            if(func != nullptr) free(func);
            func = (Float*) malloc(n * sizeof(Float));
        }

        func[0] = *f;
        // this one is not exactly correct. The constructor for func takes arguments (f, f+n), which I can't match to any
        // actual constructor, and the result isn't really what I expect. Using this one returns an image, but it is like
        // a boosted contrast
        func[1] = *f + n;

        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = Float(i) / Float(n);
        } else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }
    __both__
    int Count() const { return N-1; }
    Float SampleContinuous(Float u, Float *pdf, int *off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval(N,
                                  [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        Float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            assert(cdf[offset + 1] > cdf[offset]);
            du /= (cdf[offset + 1] - cdf[offset]);
        }
        assert(!isnan(du));

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }
    __both__
    int SampleDiscrete(Float u, Float *pdf = nullptr,
                       Float *uRemapped = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindIntervalSampling(N, cdf, u);
        // int offset = FindInterval((int)cdf.size(),
        //                           [&](int index) { return cdf[index] <= u; });
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
        if (uRemapped)
            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        if (uRemapped) assert(*uRemapped >= 0.f && *uRemapped <= 1.f);
        return offset;
    }
    Float DiscretePDF(int index) const {
        assert(index >= 0 && index < Count());
        return func[index] / (funcInt * Count());
    }

    // Distribution1D Public Data
    int N = 0;
    Float* cdf = nullptr;
    Float* func = nullptr;
    Float funcInt;
};

Point2f RejectionSampleDisk(RNG &rng);
__both__
Vector3f UniformSampleHemisphere(const Point2f &u);
__both__
Float UniformHemispherePdf();
__both__
Vector3f UniformSampleSphere(const Point2f &u);
__both__
Float UniformSpherePdf();
__both__
Vector3f UniformSampleCone(const Point2f &u, Float thetamax);
__both__
Vector3f UniformSampleCone(const Point2f &u, Float thetamax, const Vector3f &x,
                           const Vector3f &y, const Vector3f &z);
__both__
Float UniformConePdf(Float thetamax);
__both__
Point2f UniformSampleDisk(const Point2f &u);
__both__
Point2f ConcentricSampleDisk(const Point2f &u);
__both__
Point2f UniformSampleTriangle(const Point2f &u);
class Distribution2D {
  public:
    // Distribution2D Public Methods
    Distribution2D(const Float *data, int nu, int nv);
    Point2f SampleContinuous(const Point2f &u, Float *pdf) const {
        Float pdfs[2];
        int v;
        Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
        Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
        return Point2f(d0, d1);
    }
    __both__
    Float Pdf(const Point2f &p) const {
        int iu = Clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                       pConditionalV[0]->Count() - 1);
        int iv =
            Clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
    }

  private:
    // Distribution2D Private Data
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
    std::unique_ptr<Distribution1D> pMarginal;
};

// Sampling Inline Functions
template <typename T>
void Shuffle(T *samp, int count, int nDimensions, RNG &rng) {
    for (int i = 0; i < count; ++i) {
        int other = i + rng.UniformUInt32(count - i);
        for (int j = 0; j < nDimensions; ++j)
            pbrt::gpu::Swap(samp[nDimensions * i + j], samp[nDimensions * other + j]);
    }
}
__both__
inline Vector3f CosineSampleHemisphere(const Point2f &u) {
    Point2f d = ConcentricSampleDisk(u);
    Float z = sqrt(max((Float)0, 1 - d.x * d.x - d.y * d.y));
    return Vector3f(d.x, d.y, z);
}
__both__
inline Float CosineHemispherePdf(Float cosTheta) { return cosTheta * InvPi; }
__both__
inline Float BalanceHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}
__both__
inline Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    Float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

}  // namespace gpu
}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLING_H
