
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


// core/sampler.cpp*
#include "sampler.cuh"
#include "sampling.cuh"
#include "camera.cuh"
#include "stats.cuh"
#include "paramset.cuh"
#include "rng.cuh"
#include "lowdiscrepancy.cuh"

namespace pbrt {

// Sampler Local Constants
static PBRT_CONSTEXPR int kMaxResolution = 128;

// Sampler Utility Functions
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y);
static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
    int64_t x, y;
    extendedGCD(a, n, &x, &y);
    return Mod(x, n);
}

static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}

// Sampler Method Definitions
Sampler::~Sampler() {}

Sampler::Sampler(int samplesPerPixel, const Bounds2i &sampleBounds,
                             bool sampleAtPixelCenter)
    : samplesPerPixel(samplesPerPixel), sampleAtPixelCenter(sampleAtPixelCenter) {
    // Generate random digit permutations for Halton sampler
    if (radicalInversePermutations.empty()) {
        RNG rng;
        radicalInversePermutations = ComputeRadicalInversePermutations(rng);
        radicalInversePermutations_ptr = radicalInversePermutations.data();
    }

    // Find radical inverse base scales and exponents that cover sampling area
    Vector2i res = sampleBounds.pMax - sampleBounds.pMin;
    for (int i = 0; i < 2; ++i) {
        int base = (i == 0) ? 2 : 3;
        int scale = 1, exp = 0;
        while (scale < std::min(res[i], kMaxResolution)) {
            scale *= base;
            ++exp;
        }
        baseScales[i] = scale;
        baseExponents[i] = exp;
    }

    // Compute stride in samples for visiting each pixel area
    sampleStride = baseScales[0] * baseScales[1];

    // Compute multiplicative inverses for _baseScales_
    multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
    multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
}

std::vector<uint16_t> Sampler::radicalInversePermutations;
int64_t Sampler::GetIndexForSample(int64_t sampleNum) const {
    if (currentPixel != pixelForOffset) {
        // Compute Halton sample offset for _currentPixel_
        offsetForCurrentPixel = 0;
        if (sampleStride > 1) {
            Point2i pm(Mod(currentPixel[0], kMaxResolution),
                       Mod(currentPixel[1], kMaxResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                    (i == 0)
                        ? InverseRadicalInverse<2>(pm[i], baseExponents[i])
                        : InverseRadicalInverse<3>(pm[i], baseExponents[i]);
                offsetForCurrentPixel +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            offsetForCurrentPixel %= sampleStride;
        }
        pixelForOffset = currentPixel;
    }
    return offsetForCurrentPixel + sampleNum * sampleStride;
}
__both__
Float Sampler::SampleDimension(int64_t index, int dim) const {
    if (sampleAtPixelCenter && (dim == 0 || dim == 1)) return 0.5f;
    if (dim == 0)
        return RadicalInverse(dim, index >> baseExponents[0]);
    else if (dim == 1)
        return RadicalInverse(dim, index / baseScales[1]);
    else
        return ScrambledRadicalInverse(dim, index,
                                       PermutationForDimension(dim));
}

Sampler* Sampler::Clone(int seed) {
    void* samplerAddress;
    cudaMallocManaged(&samplerAddress, sizeof(Sampler));
    // LOG(ERROR) << "\n" << cudaGetErrorString(cudaGetLastError()) << std::endl;
    return new(samplerAddress) Sampler(*this);
}

CameraSample Sampler::GetCameraSample(const Point2i &pRaster) {
    CameraSample cs;
    cs.pFilm = (Point2f)pRaster + Get2D();
    cs.time = Get1D();
    cs.pLens = Get2D();
    return cs;
}

void Sampler::Request1DArray(int n) {
    CHECK_EQ(RoundCount(n), n);
    samples1DArraySizes.push_back(n);
    sampleArray1D.push_back(std::vector<Float>(n * samplesPerPixel));
}

void Sampler::Request2DArray(int n) {
    CHECK_EQ(RoundCount(n), n);
    samples2DArraySizes.push_back(n);
    sampleArray2D.push_back(std::vector<Point2f>(n * samplesPerPixel));
}

const Float *Sampler::Get1DArray(int n) {
    if (array1DOffset == sampleArray1D.size()) return nullptr;
    CHECK_EQ(samples1DArraySizes[array1DOffset], n);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return &sampleArray1D[array1DOffset++][currentPixelSampleIndex * n];
}

const Point2f *Sampler::Get2DArray(int n) {
    if (array2DOffset == sampleArray2D.size()) return nullptr;
    CHECK_EQ(samples2DArraySizes[array2DOffset], n);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
}

void Sampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    currentPixel = p;
    currentPixelSampleIndex = 0;
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(0);
    // Compute _arrayEndDim_ for dimensions used for array samples
    arrayEndDim =
        arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();

    // Compute 1D array samples for _GlobalSampler_
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        int nSamples = samples1DArraySizes[i] * samplesPerPixel;
        for (int j = 0; j < nSamples; ++j) {
            int64_t index = GetIndexForSample(j);
            sampleArray1D[i][j] = SampleDimension(index, arrayStartDim + i);
        }
    }

    // Compute 2D array samples for _GlobalSampler_
    int dim = arrayStartDim + samples1DArraySizes.size();
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        int nSamples = samples2DArraySizes[i] * samplesPerPixel;
        for (int j = 0; j < nSamples; ++j) {
            int64_t idx = GetIndexForSample(j);
            sampleArray2D[i][j].x = SampleDimension(idx, dim);
            sampleArray2D[i][j].y = SampleDimension(idx, dim + 1);
        }
        dim += 2;
    }
    CHECK_EQ(arrayEndDim, dim);
}

bool Sampler::StartNextSample() {
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(currentPixelSampleIndex + 1);
    array1DOffset = array2DOffset = 0;
    return ++currentPixelSampleIndex < samplesPerPixel;
}

bool Sampler::SetSampleNumber(int64_t sampleNum) {
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(sampleNum);
    array1DOffset = array2DOffset = 0;
    currentPixelSampleIndex = sampleNum;
    return currentPixelSampleIndex < samplesPerPixel;
}

// Float Sampler::Get1D() {
//     // ProfilePhase _(Prof::GetSample);
//     if (dimension >= arrayStartDim && dimension < arrayEndDim)
//         dimension = arrayEndDim;
//     return SampleDimension(intervalSampleIndex, dimension++);
// }

// Point2f GlobalSampler::Get2D() {
//     // ProfilePhase _(Prof::GetSample);
//     if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim)
//         dimension = arrayEndDim;
//     Point2f p(SampleDimension(intervalSampleIndex, dimension),
//               SampleDimension(intervalSampleIndex, dimension + 1));
//     dimension += 2;
//     return p;
// }

}  // namespace pbrt
