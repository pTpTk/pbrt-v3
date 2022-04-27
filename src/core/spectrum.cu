
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


// core/spectrum.cpp*
#include "spectrum.cuh"
#include <algorithm>

namespace pbrt {

// Spectrum Method Definitions
bool SpectrumSamplesSorted(const Float *lambda, const Float *vals, int n) {
    for (int i = 0; i < n - 1; ++i)
        if (lambda[i] > lambda[i + 1]) return false;
    return true;
}

void SortSpectrumSamples(Float *lambda, Float *vals, int n) {
    std::vector<std::pair<Float, Float>> sortVec;
    sortVec.reserve(n);
    for (int i = 0; i < n; ++i)
        sortVec.push_back(std::make_pair(lambda[i], vals[i]));
    std::sort(sortVec.begin(), sortVec.end());
    for (int i = 0; i < n; ++i) {
        lambda[i] = sortVec[i].first;
        vals[i] = sortVec[i].second;
    }
}

Float AverageSpectrumSamples(const Float *lambda, const Float *vals, int n,
                             Float lambdaStart, Float lambdaEnd) {
    for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
    CHECK_LT(lambdaStart, lambdaEnd);
    // Handle cases with out-of-bounds range or single sample only
    if (lambdaEnd <= lambda[0]) return vals[0];
    if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
    if (n == 1) return vals[0];
    Float sum = 0;
    // Add contributions of constant segments before/after samples
    if (lambdaStart < lambda[0]) sum += vals[0] * (lambda[0] - lambdaStart);
    if (lambdaEnd > lambda[n - 1])
        sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

    // Advance to first relevant wavelength segment
    int i = 0;
    while (lambdaStart > lambda[i + 1]) ++i;
    CHECK_LT(i + 1, n);

    // Loop over wavelength sample segments and add contributions
    auto interp = [lambda, vals](Float w, int i) {
        return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]), vals[i],
                    vals[i + 1]);
    };
    for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
        Float segLambdaStart = std::max(lambdaStart, lambda[i]);
        Float segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
        sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
               (segLambdaEnd - segLambdaStart);
    }
    return sum / (lambdaEnd - lambdaStart);
}
__both__
int FindIntervalSpectrum(int size, const Float *nodes, const Float x) {
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

Float InterpolateSpectrumSamples(const Float *lambda, const Float *vals, int n,
                                 Float l) {
    for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
    if (l <= lambda[0]) return vals[0];
    if (l >= lambda[n - 1]) return vals[n - 1];
    int offset = FindIntervalSpectrum(n, lambda, l);
    CHECK(l >= lambda[offset] && l <= lambda[offset + 1]);
    Float t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);
    return Lerp(t, vals[offset], vals[offset + 1]);
}

void Blackbody(const Float *lambda, int n, Float T, Float *Le) {
    if (T <= 0) {
        for (int i = 0; i < n; ++i) Le[i] = 0.f;
        return;
    }
    const Float c = 299792458;
    const Float h = 6.62606957e-34;
    const Float kb = 1.3806488e-23;
    for (int i = 0; i < n; ++i) {
        // Compute emitted radiance for blackbody at wavelength _lambda[i]_
        Float l = lambda[i] * 1e-9;
        Float lambda5 = (l * l) * (l * l) * l;
        Le[i] = (2 * h * c * c) /
                (lambda5 * (std::exp((h * c) / (l * kb * T)) - 1));
        CHECK(!std::isnan(Le[i]));
    }
}

void BlackbodyNormalized(const Float *lambda, int n, Float T, Float *Le) {
    Blackbody(lambda, n, T, Le);
    // Normalize _Le_ values based on maximum blackbody radiance
    Float lambdaMax = 2.8977721e-3 / T * 1e9;
    Float maxL;
    Blackbody(&lambdaMax, 1, T, &maxL);
    for (int i = 0; i < n; ++i) Le[i] /= maxL;
}

// Given a piecewise-linear SPD with values in vIn[] at corresponding
// wavelengths lambdaIn[], where lambdaIn is assumed to be sorted but may
// be irregularly spaced, resample the spectrum over the range of
// wavelengths [lambdaMin, lambdaMax], with a total of nOut wavelength
// samples between lambdaMin and lamdbaMax (including those at
// endpoints). The resampled spectrum values are written to vOut.
//
// In general, this is a standard sampling and reconstruction problem, with
// the complication that for any given invocation, some of the
// reconstruction points may involve upsampling the input distribution and
// others may involve downsampling. For upsampling, we just point-sample,
// and for downsampling, we apply a box filter centered around the
// destination wavelength with total width equal to the sample spacing.
void ResampleLinearSpectrum(const Float *lambdaIn, const Float *vIn, int nIn,
                            Float lambdaMin, Float lambdaMax, int nOut,
                            Float *vOut) {
    CHECK_GE(nOut, 2);
    for (int i = 0; i < nIn - 1; ++i) CHECK_GT(lambdaIn[i + 1], lambdaIn[i]);
    CHECK_LT(lambdaMin, lambdaMax);

    // Spacing between samples in the output distribution.
    Float delta = (lambdaMax - lambdaMin) / (nOut - 1);

    // We assume that the SPD is constant outside of the specified
    // wavelength range, taking on the respectively first/last SPD value
    // for out-of-range wavelengths.
    //
    // To make this convention fit cleanly into the code below, we create
    // virtual samples in the input distribution with index -1 for the
    // sample before the first valid sample and index nIn for the sample
    // after the last valid sample. In turn, can place those virtual
    // samples beyond the endpoints of the target range so that we can
    // always assume that the source range is broader than the target
    // range, which in turn lets us not worry about various boundary cases
    // below.

    // The wavelengths of the virtual samples at the endpoints are set so
    // that they are one destination sample spacing beyond the destination
    // endpoints.  (Note that this potentially means that if we swept along
    // indices from -1 to nIn, we wouldn't always see a monotonically
    // increasing set of wavelength values. However, this isn't a problem
    // since we only access these virtual samples if the destination range
    // is wider than the source range.)
    auto lambdaInClamped = [&](int index) {
        CHECK(index >= -1 && index <= nIn);
        if (index == -1) {
            CHECK_LT(lambdaMin - delta, lambdaIn[0]);
            return lambdaMin - delta;
        } else if (index == nIn) {
            CHECK_GT(lambdaMax + delta, lambdaIn[nIn - 1]);
            return lambdaMax + delta;
        }
        return lambdaIn[index];
    };

    // Due to the piecewise-constant assumption, the SPD values outside the
    // specified range are given by the valid endpoints.
    auto vInClamped = [&](int index) {
        CHECK(index >= -1 && index <= nIn);
        return vIn[Clamp(index, 0, nIn - 1)];
    };

    // Helper that upsamples ors downsample the given SPD at the given
    // wavelength lambda.
    auto resample = [&](Float lambda) -> Float {
        // Handle the edge cases first so that we don't need to worry about
        // them in the following.
        //
        // First, if the entire filtering range for the destination is
        // outside of the range of valid samples, we can just return the
        // endpoint value.
        if (lambda + delta / 2 <= lambdaIn[0]) return vIn[0];
        if (lambda - delta / 2 >= lambdaIn[nIn - 1]) return vIn[nIn - 1];
        // Second, if there's only one sample, then the SPD has the same
        // value everywhere, and we're done.
        if (nIn == 1) return vIn[0];

        // Otherwise, find indices into the input SPD that bracket the
        // wavelength range [lambda-delta, lambda+delta]. Note that this is
        // a 2x wider range than we will actually filter over in the end.
        int start, end;
        if (lambda - delta < lambdaIn[0])
            // Virtual sample at the start, as described above.
            start = -1;
        else {
            start = FindInterval(
                nIn, [&](int i) { return lambdaIn[i] <= lambda - delta; });
            CHECK(start >= 0 && start < nIn);
        }

        if (lambda + delta > lambdaIn[nIn - 1])
            // Virtual sample at the end, as described above.
            end = nIn;
        else {
            // Linear search from the starting point. (Presumably more
            // efficient than a binary search from scratch, or doesn't
            // matter either way.)
            end = start > 0 ? start : 0;
            while (end < nIn && lambda + delta > lambdaIn[end]) ++end;
        }

        if (end - start == 2 && lambdaInClamped(start) <= lambda - delta &&
            lambdaIn[start + 1] == lambda &&
            lambdaInClamped(end) >= lambda + delta) {
            // Downsampling: special case where the input and output
            // wavelengths line up perfectly, so just return the
            // corresponding point sample at lambda.
            return vIn[start + 1];
        } else if (end - start == 1) {
            // Downsampling: evaluate the piecewise-linear function at
            // lambda.
            Float t = (lambda - lambdaInClamped(start)) /
                      (lambdaInClamped(end) - lambdaInClamped(start));
            CHECK(t >= 0 && t <= 1);
            return Lerp(t, vInClamped(start), vInClamped(end));
        } else {
            // Upsampling: use a box filter and average all values in the
            // input spectrum from lambda +/- delta / 2.
            return AverageSpectrumSamples(
                lambdaIn, vIn, nIn, lambda - delta / 2, lambda + delta / 2);
        }
    };

    // For each destination sample, compute the wavelength lambda for the
    // sample and then resample the source SPD distribution at that point.
    for (int outOffset = 0; outOffset < nOut; ++outOffset) {
        // TODO: Currently, resample() does a binary search each time,
        // even though we could do a single sweep across the input array,
        // since we're resampling it at a regular and increasing set of
        // lambdas. It would be nice to polish that up.
        Float lambda =
            Lerp(Float(outOffset) / (nOut - 1), lambdaMin, lambdaMax);
        vOut[outOffset] = resample(lambda);
    }
}

}  // namespace pbrt
