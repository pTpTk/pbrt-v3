
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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.cuh"
#include "geometry.cuh"
#include "rng.cuh"
#include <inttypes.h>
#include "lowdiscrepancy.cuh"

namespace pbrt {

// Sampler Declarations
class Sampler {
  public:
    // Sampler Interface
    ~Sampler();
    Sampler(int64_t samplesPerPixel);
    Sampler(int nsamp, const Bounds2i &sampleBounds,
            bool sampleAtCenter = false);
    void StartPixel(const Point2i &p);
    __both__
    Float Get1D();
    __both__
    Point2f Get2D();
    CameraSample GetCameraSample(const Point2i &pRaster);
    void Request1DArray(int n);
    void Request2DArray(int n);
    int RoundCount(int n) const { return n; }
    const Float *Get1DArray(int n);
    const Point2f *Get2DArray(int n);
    bool StartNextSample();
    Sampler* Clone(int seed);
    bool SetSampleNumber(int64_t sampleNum);
    std::string StateString() const {
      return StringPrintf("(%d,%d), sample %" PRId64, currentPixel.x,
                          currentPixel.y, currentPixelSampleIndex);
    }
    int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }

    int64_t GetIndexForSample(int64_t sampleNum) const;
    __both__
    Float SampleDimension(int64_t index, int dimension) const;

    // Sampler Public Data
    const int64_t samplesPerPixel;

  protected:
    // Sampler Protected Data
    Point2i currentPixel;
    int64_t currentPixelSampleIndex;
    std::vector<int> samples1DArraySizes, samples2DArraySizes;
    std::vector<std::vector<Float>> sampleArray1D;
    std::vector<std::vector<Point2f>> sampleArray2D;

  private:
    // Sampler Private Data
    size_t array1DOffset, array2DOffset;
    
    // GlobalSampler Private Data
    int dimension;
    int64_t intervalSampleIndex;
    static const int arrayStartDim = 5;
    int arrayEndDim;

    // HaltonSampler Private Data
    static std::vector<uint16_t> radicalInversePermutations;
    uint16_t* radicalInversePermutations_ptr;
    Point2i baseScales, baseExponents;
    int sampleStride;
    int multInverse[2];
    mutable Point2i pixelForOffset = Point2i(numeric_limits<int>::max(),
                                             numeric_limits<int>::max());
    mutable int64_t offsetForCurrentPixel;
    // Added after book publication: force all image samples to be at the
    // center of the pixel area.
    bool sampleAtPixelCenter;

    // HaltonSampler Private Methods
    __both__
    const uint16_t *PermutationForDimension(int dim) const {
        return &radicalInversePermutations_ptr[PrimeSums[dim]];
    }
};

// class GlobalSampler : public Sampler {
//   public:
//     // GlobalSampler Public Methods
//     bool StartNextSample();
//     void StartPixel(const Point2i &);
//     bool SetSampleNumber(int64_t sampleNum);
//     __both__
//     Float Get1D();
//     __both__
//     Point2f Get2D();
//     GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}
//     virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
//     __both__
//     virtual Float SampleDimension(int64_t index, int dimension) const = 0;

//   private:
//     // GlobalSampler Private Data
//     int dimension;
//     int64_t intervalSampleIndex;
//     static const int arrayStartDim = 5;
//     int arrayEndDim;
// };

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLER_H
