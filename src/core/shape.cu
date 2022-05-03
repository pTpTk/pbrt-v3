
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


// core/shape.cpp*
#include "shape.cuh"
#include "stats.cuh"
#include "lowdiscrepancy.cuh"
#include "sampling.cuh"
#include "paramset.cuh"
#include "efloat.cuh"

namespace pbrt {

// Shape Method Definitions
Shape::~Shape() {}

STAT_COUNTER("Scene/Shapes created", nShapesCreated);
// Shape::Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
//              bool reverseOrientation)
//     : ObjectToWorld(ObjectToWorld),
//       WorldToObject(WorldToObject),
//       reverseOrientation(reverseOrientation),
//       transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
//     ++nShapesCreated;
// }

Shape::Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
           bool reverseOrientation, Float radius, Float zMin, Float zMax,
           Float phiMax)
        : ObjectToWorld(ObjectToWorld),
          WorldToObject(WorldToObject),
          reverseOrientation(reverseOrientation),
          transformSwapsHandedness(ObjectToWorld->SwapsHandedness()),
          radius(radius),
          zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
          zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
          thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
          thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
          phiMax(Radians(Clamp(phiMax, 0, 360))) {
              ++nShapesCreated;
          }

// __both__
// Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
// __both__
// Interaction Shape::Sample(const Interaction &ref, const Point2f &u,
//                           Float *pdf) const {
//     Interaction intr = Sample(u, pdf);
//     Vector3f wi = intr.p - ref.p;
//     if (wi.LengthSquared() == 0)
//         *pdf = 0;
//     else {
//         wi = Normalize(wi);
//         // Convert from area measure, as returned by the Sample() call
//         // above, to solid angle measure.
//         *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
//         if (isinf(*pdf)) *pdf = 0.f;
//     }
//     return intr;
// }
// __both__
// Float Shape::Pdf(const Interaction &ref, const Vector3f &wi) const {
//     // Intersect sample ray with area light geometry
//     Ray ray = ref.SpawnRay(wi);
//     Float tHit;
//     SurfaceInteraction isectLight;
//     // Ignore any alpha textures used for trimming the shape when performing
//     // this intersection. Hack for the "San Miguel" scene, where this is used
//     // to make an invisible area light.
//     if (!Intersect(ray, &tHit, &isectLight, false)) return 0;

//     // Convert light sample weight to solid angle measure
//     Float pdf = DistanceSquared(ref.p, isectLight.p) /
//                 (AbsDot(isectLight.n, -wi) * Area());
//     if (isinf(pdf)) pdf = 0.f;
//     return pdf;
// }

// Float Shape::SolidAngle(const Point3f &p, int nSamples) const {
//     Interaction ref(p, Normal3f(), Vector3f(), Vector3f(0, 0, 1), 0,
//                     MediumInterface{});
//     double solidAngle = 0;
//     for (int i = 0; i < nSamples; ++i) {
//         Point2f u{RadicalInverse(0, i), RadicalInverse(1, i)};
//         Float pdf;
//         Interaction pShape = Sample(ref, u, &pdf);
//         if (pdf > 0 && !IntersectP(Ray(p, pShape.p - p, .999f))) {
//             solidAngle += 1 / pdf;
//         }
//     }
//     return solidAngle / nSamples;
// }

__both__
Bounds3f Shape::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, zMin),
                    Point3f(radius, radius, zMax));
}
__both__
bool Shape::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect) const {
    // ProfilePhase p(Prof::ShapeIntersect);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute sphere hit position and $\phi$
    pHit = ray((Float)tShapeHit);

    // Refine sphere intersection point
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > ray.tMax) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    // Find parametric representation of sphere hit
    Float u = phi / phiMax;
    Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    Float zRadius = pbrt::math::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    Float invZRadius = 1 / zRadius;
    Float cosPhi = pHit.x * invZRadius;
    Float sinPhi = pHit.y * invZRadius;
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv =
        (thetaMax - thetaMin) *
        Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

    // Compute sphere $\dndu$ and $\dndv$
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv =
        (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                      Vector3f(pHit.x, pHit.y, pHit.z);

    // Compute coefficients for fundamental forms
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    Float invEGF2 = 1 / (E * G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
                             (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
                             (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for sphere intersection
    Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

    // Initialize _SurfaceInteraction_ from parametric information
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;
    return true;
}
__both__
bool Shape::IntersectP(const Ray &r) const {
    // ProfilePhase p(Prof::ShapeIntersectP);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute sphere hit position and $\phi$
    pHit = ray((Float)tShapeHit);

    // Refine sphere intersection point
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > ray.tMax) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }
    return true;
}
__both__
Float Shape::Area() const { return phiMax * radius * (zMax - zMin); }
__both__
Interaction Shape::Sample(const Point2f &u, Float *pdf) const {
    Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
    if (reverseOrientation) it.n *= -1;
    // Reproject _pObj_ to sphere surface and compute _pObjError_
    pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
    Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
    it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
    *pdf = 1 / Area();
    return it;
}
__both__
Interaction Shape::Sample(const Interaction &ref, const Point2f &u,
                           Float *pdf) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

    // Sample uniformly on sphere if $\pt{}$ is inside it
    Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
    if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
        Interaction intr = Sample(u, pdf);
        Vector3f wi = intr.p - ref.p;
        if (wi.LengthSquared() == 0)
            *pdf = 0;
        else {
            // Convert from area measure returned by Sample() call above to
            // solid angle measure.
            wi = Normalize(wi);
            *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
        }
        if (isinf(*pdf)) *pdf = 0.f;
        return intr;
    }

    // Sample sphere uniformly inside subtended cone

    // Compute coordinate system for sphere sampling
    Float dc = Distance(ref.p, pCenter);
    Float invDc = 1 / dc;
    Vector3f wc = (pCenter - ref.p) * invDc;
    Vector3f wcX, wcY;
    CoordinateSystem(wc, &wcX, &wcY);

    // Compute $\theta$ and $\phi$ values for sample in cone
    Float sinThetaMax = radius * invDc;
    Float sinThetaMax2 = sinThetaMax * sinThetaMax;
    Float invSinThetaMax = 1 / sinThetaMax;
    Float cosThetaMax = pbrt::math::sqrt(max((Float)0.f, 1 - sinThetaMax2));

    Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
    Float sinTheta2 = 1 - cosTheta * cosTheta;

    if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
        /* Fall back to a Taylor series expansion for small angles, where
           the standard approach suffers from severe cancellation errors */
        sinTheta2 = sinThetaMax2 * u[0];
        cosTheta = pbrt::math::sqrt(1 - sinTheta2);
    }

    // Compute angle $\alpha$ from center of sphere to sampled point on surface
    Float cosAlpha = sinTheta2 * invSinThetaMax +
        cosTheta * pbrt::math::sqrt(max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
    Float sinAlpha = pbrt::math::sqrt(max((Float)0.f, 1.f - cosAlpha*cosAlpha));
    Float phi = u[1] * 2 * Pi;

    // Compute surface normal and sampled point on sphere
    Vector3f nWorld =
        SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
    Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

    // Return _Interaction_ for sampled point on sphere
    Interaction it;
    it.p = pWorld;
    it.pError = gamma(5) * Abs((Vector3f)pWorld);
    it.n = Normal3f(nWorld);
    if (reverseOrientation) it.n *= -1;

    // Uniform cone PDF.
    *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

    return it;
}
__both__
Float Shape::Pdf(const Interaction &ref, const Vector3f &wi) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
    // Return uniform PDF if point is inside sphere
    Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
    if (DistanceSquared(pOrigin, pCenter) <= radius * radius){
        // Intersect sample ray with area light geometry
        Ray ray = ref.SpawnRay(wi);
        Float tHit;
        SurfaceInteraction isectLight;
        // Ignore any alpha textures used for trimming the shape when performing
        // this intersection. Hack for the "San Miguel" scene, where this is used
        // to make an invisible area light.
        if (!Intersect(ray, &tHit, &isectLight)) return 0;

        // Convert light sample weight to solid angle measure
        Float pdf = DistanceSquared(ref.p, isectLight.p) /
                    (AbsDot(isectLight.n, -wi) * Area());
        if (isinf(pdf)) pdf = 0.f;
        return pdf;
    }

    // Compute general sphere PDF
    Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
    Float cosThetaMax = pbrt::math::sqrt(max((Float)0, 1 - sinThetaMax2));
    return UniformConePdf(cosThetaMax);
}

Float Shape::SolidAngle(const Point3f &p, int nSamples) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
    if (DistanceSquared(p, pCenter) <= radius * radius)
        return 4 * Pi;
    Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
    Float cosTheta = pbrt::math::sqrt(std::max((Float)0, 1 - sinTheta2));
    return (2 * Pi * (1 - cosTheta));
}

}  // namespace pbrt
