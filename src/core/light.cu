
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


// core/light.cpp*
#include "light.cuh"
#include "scene.cuh"
#include "sampling.cuh"
#include "stats.cuh"
#include "paramset.cuh"

namespace pbrt {

STAT_COUNTER("Scene/Lights", numLights);
STAT_COUNTER("Scene/AreaLights", numAreaLights);

// Light Method Definitions
// Light::Light(int flags, const Transform &LightToWorld,
//              const MediumInterface &mediumInterface, int nSamples)
//     : flags(flags),
//       nSamples(std::max(1, nSamples)),
//       mediumInterface(mediumInterface),
//       LightToWorld(LightToWorld),
//       WorldToLight(Inverse(LightToWorld)) {
//     ++numLights;
// }

Light::~Light() {}
// __both__
// bool VisibilityTester::Unoccluded(const Scene &scene) const {
//     return !scene.IntersectP(p0.SpawnRayTo(p1));
// }
// __both__
// Spectrum Light::Le(const RayDifferential &ray) const { return Spectrum(0.f); }

// Light::Light(const Transform &LightToWorld, const MediumInterface &medium,
//                      int nSamples)
//     : flags((int)LightFlags::Area),
//       nSamples(std::max(1, nSamples)),
//       mediumInterface(medium),
//       LightToWorld(LightToWorld),
//       WorldToLight(Inverse(LightToWorld)) {
//     ++numAreaLights;
// }

Light::Light(const Transform &LightToWorld,
             const MediumInterface &mediumInterface,
             const Spectrum &Lemit, int nSamples,
             Shape* const shape,
             bool twoSided)
    : flags((int)LightFlags::Area),
      nSamples(std::max(1, nSamples)),
      mediumInterface(mediumInterface),
      LightToWorld(LightToWorld),
      WorldToLight(Inverse(LightToWorld)),
      Lemit(Lemit),
      shape(shape),
      twoSided(twoSided),
      area(shape->Area()) {
        ++numLights;
        ++numAreaLights;
    // Warn if light has transformation with non-uniform scale, though not
    // for Triangles, since this doesn't matter for them.
    if (WorldToLight.HasScale() &&
        shape == nullptr)
        Warning(
            "Scaling detected in world to light transformation! "
            "The system has numerous assumptions, implicit and explicit, "
            "that this transform will have no scale factors in it. "
            "Proceed at your own risk; your image may have errors.");
}

Spectrum Light::Power() const {
    return (twoSided ? 2 : 1) * Lemit * area * Pi;
}
__both__
Spectrum Light::Sample_Li(const Interaction &ref, const Point2f &u,
                                     Vector3f *wi, Float *pdf,
                                     VisibilityTester *vis) const {
    // ProfilePhase _(Prof::LightSample);
    Interaction pShape = shape->Sample(ref, u, pdf);
    pShape.mediumInterface = mediumInterface;
    if (*pdf == 0 || (pShape.p - ref.p).LengthSquared() == 0) {
        *pdf = 0;
        return 0.f;
    }
    *wi = Normalize(pShape.p - ref.p);
    *vis = VisibilityTester(ref, pShape);
    return L(pShape, -*wi);
}
__both__
Float Light::Pdf_Li(const Interaction &ref,
                               const Vector3f &wi) const {
    // ProfilePhase _(Prof::LightPdf);
    return shape->Pdf(ref, wi);
}

Spectrum Light::Sample_Le(const Point2f &u1, const Point2f &u2,
                                     Float time, Ray *ray, Normal3f *nLight,
                                     Float *pdfPos, Float *pdfDir) const {
    ProfilePhase _(Prof::LightSample);
    // Sample a point on the area light's _Shape_, _pShape_
    Interaction pShape = shape->Sample(u1, pdfPos);
    pShape.mediumInterface = mediumInterface;
    *nLight = pShape.n;

    // Sample a cosine-weighted outgoing direction _w_ for area light
    Vector3f w;
    if (twoSided) {
        Point2f u = u2;
        // Choose a side to sample and then remap u[0] to [0,1] before
        // applying cosine-weighted hemisphere sampling for the chosen side.
        if (u[0] < .5) {
            u[0] = std::min(u[0] * 2, OneMinusEpsilon);
            w = CosineSampleHemisphere(u);
        } else {
            u[0] = std::min((u[0] - .5f) * 2, OneMinusEpsilon);
            w = CosineSampleHemisphere(u);
            w.z *= -1;
        }
        *pdfDir = 0.5f * CosineHemispherePdf(std::abs(w.z));
    } else {
        w = CosineSampleHemisphere(u2);
        *pdfDir = CosineHemispherePdf(w.z);
    }

    Vector3f v1, v2, n(pShape.n);
    CoordinateSystem(n, &v1, &v2);
    w = w.x * v1 + w.y * v2 + w.z * n;
    *ray = pShape.SpawnRay(w);
    return L(pShape, w);
}

void Light::Pdf_Le(const Ray &ray, const Normal3f &n, Float *pdfPos,
                              Float *pdfDir) const {
    ProfilePhase _(Prof::LightPdf);
    Interaction it(ray.o, n, Vector3f(), Vector3f(n), ray.time,
                   mediumInterface);
    *pdfPos = shape->Pdf(it);
    *pdfDir = twoSided ? (.5 * CosineHemispherePdf(AbsDot(n, ray.d)))
                       : CosineHemispherePdf(Dot(n, ray.d));
}

}  // namespace pbrt
