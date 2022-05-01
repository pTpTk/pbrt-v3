
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

#ifndef PBRT_SAMPLERS_HALTON_H
#define PBRT_SAMPLERS_HALTON_H

// samplers/halton.h*
#include "sampler.cuh"
#include "lowdiscrepancy.cuh"

namespace pbrt {

// HaltonSampler Declarations
// class HaltonSampler : public GlobalSampler {
//   public:
//     // HaltonSampler Public Methods
//     HaltonSampler(int nsamp, const Bounds2i &sampleBounds,
//                   bool sampleAtCenter = false);
//     int64_t GetIndexForSample(int64_t sampleNum) const;
//     __both__
//     Float SampleDimension(int64_t index, int dimension) const;
//     Sampler* Clone(int seed);

//   private:
//     // HaltonSampler Private Data
//     static std::vector<uint16_t> radicalInversePermutations;
//     uint16_t* radicalInversePermutations_ptr;
//     Point2i baseScales, baseExponents;
//     int sampleStride;
//     int multInverse[2];
//     mutable Point2i pixelForOffset = Point2i(numeric_limits<int>::max(),
//                                              numeric_limits<int>::max());
//     mutable int64_t offsetForCurrentPixel;
//     // Added after book publication: force all image samples to be at the
//     // center of the pixel area.
//     bool sampleAtPixelCenter;

//     // HaltonSampler Private Methods
//     __both__
//     const uint16_t *PermutationForDimension(int dim) const {
//         return &radicalInversePermutations_ptr[PrimeSums[dim]];
//     }
// };

Sampler *CreateHaltonSampler(const ParamSet &params,
                                   const Bounds2i &sampleBounds);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_HALTON_H
