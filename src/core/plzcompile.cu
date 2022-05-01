// core/api.cpp*
#include "api.cuh"
#include "parallel.cuh"
#include "paramset.cuh"
#include "spectrum.cuh"
#include "scene.cuh"
#include "film.cuh"
#include "medium.cuh"
#include "stats.cuh"

// API Additional Headers
#include "accelerators/bvh.cuh"
#include "cameras/perspective.cuh"
#include "filters/box.cuh"
#include "integrators/path.cuh"
#include "lights/diffuse.cuh"
#include "materials/matte.cuh"
#include "samplers/halton.cuh"
#include "shapes/sphere.cuh"
#include "progressreporter.cuh"

namespace pbrt{

// integrators/path
__device__
Spectrum PathIntegrator::Li(const RayDifferential &r, const Scene &scene,
                            Sampler &sampler, MemoryArena &arena,
                            int depth) const {
    // ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
    int bounces;
    // Added after book publication: etaScale tracks the accumulated effect
    // of radiance scaling due to rays passing through refractive
    // boundaries (see the derivation on p. 527 of the third edition). We
    // track this value in order to remove it from beta when we apply
    // Russian roulette; this is worthwhile, since it lets us sometimes
    // avoid terminating refracted rays that are about to be refracted back
    // out of a medium and thus have their beta value increased.
    Float etaScale = 1;

    //for (bounces = 0;; ++bounces) {
    for (bounces = 0; bounces < 1; ++bounces) {
        // Find next path vertex and accumulate contribution

        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect;
        printf("isect[%p]\n", &isect);
        printf("scene[%p]\n", &scene);
        bool foundIntersection = scene.Intersect(ray, &isect);

        // // Possibly add emitted light at intersection
        // if (bounces == 0 || specularBounce) {
        //     // Add emitted light at path vertex or from the environment
        //     if (foundIntersection) {
        //         L += beta * isect.Le(-ray.d);
        //         // VLOG(2) << "Added Le -> L = " << L;
        //     } else {
        //         for (int i = 0; i < scene.infiniteLights_size; ++i) {
        //             const auto &light = scene.infiniteLights[i];
        //             L += beta * light->Le(ray);
        //         }
        //         // VLOG(2) << "Added infinite area lights -> L = " << L;
        //     }
        // }

        // // Terminate path if ray escaped or _maxDepth_ was reached
        // if (!foundIntersection || bounces >= maxDepth) break;

        // // Compute scattering functions and skip over medium boundaries
        // isect.ComputeScatteringFunctions(ray, arena, true);
        // if (!isect.bsdf) {
        //     // VLOG(2) << "Skipping intersection due to null bsdf";
        //     ray = isect.SpawnRay(ray.d);
        //     bounces--;
        //     continue;
        // }

        // const Distribution1D *distrib = lightDistribution->Lookup(isect.p);

        // // Sample illumination from lights to find path contribution.
        // // (But skip this for perfectly specular BSDFs.)
        // if (isect.bsdf->NumComponents(BxDFType(BSDF_ALL & ~BSDF_SPECULAR)) >
        //     0) {
        //     // ++totalPaths;
        //     Spectrum Ld = beta * UniformSampleOneLight(isect, scene, arena,
        //                                                sampler, false, distrib);
        //     // VLOG(2) << "Sampled direct lighting Ld = " << Ld;
        //     // if (Ld.IsBlack()) ++zeroRadiancePaths;
        //     // CHECK_GE(Ld.y(), 0.f);
        //     L += Ld;
        // }

        // // Sample BSDF to get new path direction
        // Vector3f wo = -ray.d, wi;
        // Float pdf;
        // BxDFType flags;
        // Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
        //                                   BSDF_ALL, &flags);
        // // VLOG(2) << "Sampled BSDF, f = " << f << ", pdf = " << pdf;
        // if (f.IsBlack() || pdf == 0.f) break;
        // beta *= f * AbsDot(wi, isect.shading.n) / pdf;
        // // VLOG(2) << "Updated beta = " << beta;
        // // CHECK_GE(beta.y(), 0.f);
        // // DCHECK(!std::isinf(beta.y()));
        // specularBounce = (flags & BSDF_SPECULAR) != 0;
        // if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
        //     Float eta = isect.bsdf->eta;
        //     // Update the term that tracks radiance scaling for refraction
        //     // depending on whether the ray is entering or leaving the
        //     // medium.
        //     etaScale *= (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
        // }
        // ray = isect.SpawnRay(wi);

        // // Possibly terminate the path with Russian roulette.
        // // Factor out radiance scaling due to refraction in rrBeta.
        // Spectrum rrBeta = beta * etaScale;
        // if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
        //     Float q = max((Float).05, 1 - rrBeta.MaxComponentValue());
        //     if (sampler.Get1D() < q) break;
        //     beta /= 1 - q;
        //     // DCHECK(!std::isinf(beta.y()));
        // }
    }
    //ReportValue(pathLength, bounces);
    return L;
}

__both__
bool Scene::IntersectP(const Ray &ray) const {
    // ++nShadowTests;
    assert(ray.d != Vector3f(0,0,0));
    return aggregate->IntersectP(ray);
}

__both__
Vector3f UniformSampleHemisphere(const Point2f &u) {
    Float z = u[0];
    Float r = pbrt::math::sqrt(max((Float)0, (Float)1. - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}
__both__
Float UniformHemispherePdf() { return Inv2Pi; }

__both__
Point2f ConcentricSampleDisk(const Point2f &u) {
    // Map uniform random numbers to $[-1,1]^2$
    Point2f uOffset = 2.f * u - Vector2f(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return Point2f(0, 0);

    // Apply concentric mapping to point
    Float theta, r;
    if (pbrt::math::abs(uOffset.x) > pbrt::math::abs(uOffset.y)) {
        r = uOffset.x;
        theta = PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
    }
    return r * Point2f(std::cos(theta), std::sin(theta));
}

__both__
void SurfaceInteraction::ComputeDifferentials(
    const RayDifferential &ray) const {
    if (ray.hasDifferentials) {
        // Estimate screen space change in $\pt{}$ and $(u,v)$

        // Compute auxiliary intersection points with plane
        Float d = Dot(n, Vector3f(p.x, p.y, p.z));
        Float tx =
            -(Dot(n, Vector3f(ray.rxOrigin)) - d) / Dot(n, ray.rxDirection);
        if (isinf(tx) || pbrt::math::isnan(tx)) goto fail;
        Point3f px = ray.rxOrigin + tx * ray.rxDirection;
        Float ty =
            -(Dot(n, Vector3f(ray.ryOrigin)) - d) / Dot(n, ray.ryDirection);
        if (isinf(ty) || pbrt::math::isnan(ty)) goto fail;
        Point3f py = ray.ryOrigin + ty * ray.ryDirection;
        dpdx = px - p;
        dpdy = py - p;

        // Compute $(u,v)$ offsets at auxiliary points

        // Choose two dimensions to use for ray offset computation
        int dim[2];
        if (pbrt::math::abs(n.x) > pbrt::math::abs(n.y) && pbrt::math::abs(n.x) > pbrt::math::abs(n.z)) {
            dim[0] = 1;
            dim[1] = 2;
        } else if (pbrt::math::abs(n.y) > pbrt::math::abs(n.z)) {
            dim[0] = 0;
            dim[1] = 2;
        } else {
            dim[0] = 0;
            dim[1] = 1;
        }

        // Initialize _A_, _Bx_, and _By_ matrices for offset computation
        Float A[2][2] = {{dpdu[dim[0]], dpdv[dim[0]]},
                         {dpdu[dim[1]], dpdv[dim[1]]}};
        Float Bx[2] = {px[dim[0]] - p[dim[0]], px[dim[1]] - p[dim[1]]};
        Float By[2] = {py[dim[0]] - p[dim[0]], py[dim[1]] - p[dim[1]]};
        if (!SolveLinearSystem2x2(A, Bx, &dudx, &dvdx)) dudx = dvdx = 0;
        if (!SolveLinearSystem2x2(A, By, &dudy, &dvdy)) dudy = dvdy = 0;
    } else {
    fail:
        dudx = dvdx = 0;
        dudy = dvdy = 0;
        dpdx = dpdy = Vector3f(0, 0, 0);
    }
}

__both__
void SurfaceInteraction::SetShadingGeometry(const Vector3f &dpdus,
                                            const Vector3f &dpdvs,
                                            const Normal3f &dndus,
                                            const Normal3f &dndvs,
                                            bool orientationIsAuthoritative) {
    // Compute _shading.n_ for _SurfaceInteraction_
    shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
    if (orientationIsAuthoritative)
        n = Faceforward(n, shading.n);
    else
        shading.n = Faceforward(shading.n, n);

    // Initialize _shading_ partial derivative values
    shading.dpdu = dpdus;
    shading.dpdv = dpdvs;
    shading.dndu = dndus;
    shading.dndv = dndvs;
}

Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

__both__
Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }

__both__
Float Shape::Pdf(const Interaction &ref, const Vector3f &wi) const {
    // Intersect sample ray with area light geometry
    Ray ray = ref.SpawnRay(wi);
    Float tHit;
    SurfaceInteraction isectLight;
    // Ignore any alpha textures used for trimming the shape when performing
    // this intersection. Hack for the "San Miguel" scene, where this is used
    // to make an invisible area light.
    if (!Intersect(ray, &tHit, &isectLight, false)) return 0;

    // Convert light sample weight to solid angle measure
    Float pdf = DistanceSquared(ref.p, isectLight.p) /
                (AbsDot(isectLight.n, -wi) * Area());
    if (isinf(pdf)) pdf = 0.f;
    return pdf;
}

__both__
Float UniformConePdf(Float cosThetaMax) {
    return 1 / (2 * Pi * (1 - cosThetaMax));
}

// primitive
__both__
const AreaLight *Aggregate::GetAreaLight() const {
    return nullptr;
}
__both__
const Material *Aggregate::GetMaterial() const {
    return nullptr;
}

__both__
void Aggregate::ComputeScatteringFunctions(SurfaceInteraction *isect,
                                           MemoryArena &arena,
                                           TransportMode mode,
                                           bool allowMultipleLobes) const {
}

__both__
Bounds3f Transform::operator()(const Bounds3f &b) const {
    const Transform &M = *this;
    Bounds3f ret(M(Point3f(b.pMin.x, b.pMin.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMin.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMax.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMin.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMax.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMax.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMin.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMax.y, b.pMax.z)));
    return ret;
}

__both__
SurfaceInteraction Transform::operator()(const SurfaceInteraction &si) const {
    SurfaceInteraction ret;
    // Transform _p_ and _pError_ in _SurfaceInteraction_
    ret.p = (*this)(si.p, si.pError, &ret.pError);

    // Transform remaining members of _SurfaceInteraction_
    const Transform &t = *this;
    ret.n = Normalize(t(si.n));
    ret.wo = Normalize(t(si.wo));
    ret.time = si.time;
    ret.mediumInterface = si.mediumInterface;
    ret.uv = si.uv;
    ret.shape = si.shape;
    ret.dpdu = t(si.dpdu);
    ret.dpdv = t(si.dpdv);
    ret.dndu = t(si.dndu);
    ret.dndv = t(si.dndv);
    ret.shading.n = Normalize(t(si.shading.n));
    ret.shading.dpdu = t(si.shading.dpdu);
    ret.shading.dpdv = t(si.shading.dpdv);
    ret.shading.dndu = t(si.shading.dndu);
    ret.shading.dndv = t(si.shading.dndv);
    ret.dudx = si.dudx;
    ret.dvdx = si.dvdx;
    ret.dudy = si.dudy;
    ret.dvdy = si.dvdy;
    ret.dpdx = t(si.dpdx);
    ret.dpdy = t(si.dpdy);
    ret.bsdf = si.bsdf;
    ret.primitive = si.primitive;
    //    ret.n = Faceforward(ret.n, ret.shading.n);
    ret.shading.n = Faceforward(ret.shading.n, ret.n);
    ret.faceIndex = si.faceIndex;
    return ret;
}

__both__
SurfaceInteraction::SurfaceInteraction(
    const Point3f &p, const Vector3f &pError, const Point2f &uv,
    const Vector3f &wo, const Vector3f &dpdu, const Vector3f &dpdv,
    const Normal3f &dndu, const Normal3f &dndv, Float time, const Shape *shape,
    int faceIndex)
    : Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), pError, wo, time,
                  nullptr),
      uv(uv),
      dpdu(dpdu),
      dpdv(dpdv),
      dndu(dndu),
      dndv(dndv),
      shape(shape),
      faceIndex(faceIndex) {
    // Initialize shading geometry from true geometry
    shading.n = n;
    shading.dpdu = dpdu;
    shading.dpdv = dpdv;
    shading.dndu = dndu;
    shading.dndv = dndv;

    // Adjust normal based on orientation and handedness
    if (shape &&
        (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
        n *= -1;
        shading.n *= -1;
    }
}

__both__
Vector3f UniformSampleSphere(const Point2f &u) {
    Float z = 1 - 2 * u[0];
    Float r = pbrt::math::sqrt(max((Float)0, (Float)1 - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

// integrator
__both__
Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, bool handleMedia, bool specular) {
    BxDFType bsdfFlags =
        specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
    Spectrum Ld(0.f);
    // Sample light source with multiple importance sampling
    Vector3f wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
    if (lightPdf > 0 && !Li.IsBlack()) {
        // Compute BSDF or phase function's value for light sample
        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            // Evaluate BSDF for light sampling strategy
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
                AbsDot(wi, isect.shading.n);
            scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
        }
        if (!f.IsBlack()) {
            // Compute effect of visibility for light source sample
              if (!visibility.Unoccluded(scene)) {
                Li = Spectrum(0.f);
            }

            // Add light's contribution to reflected radiance
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags))
                    Ld += f * Li / lightPdf;
                else {
                    Float weight =
                        PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                    Ld += f * Li * weight / lightPdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!IsDeltaLight(light.flags)) {
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction()) {
            // Sample scattered direction for surface interactions
            BxDFType sampledType;
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
                                     bsdfFlags, &sampledType);
            f *= AbsDot(wi, isect.shading.n);
            sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
        }
        if (!f.IsBlack() && scatteringPdf > 0) {
            // Account for light contributions along sampled direction _wi_
            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);
                if (lightPdf == 0) return Ld;
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }

            // Find intersection and compute transmittance
            SurfaceInteraction lightIsect;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction =
                handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
                            : scene.Intersect(ray, &lightIsect);

            // Add light contribution from material sampling
            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIsect.primitive->GetAreaLight() == &light)
                    Li = lightIsect.Le(-wi);
            } else
                Li = light.Le(ray);
            if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
        }
    }
    return Ld;
}

__both__
Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
                               MemoryArena &arena, Sampler &sampler,
                               bool handleMedia, const Distribution1D *lightDistrib) {
    // ProfilePhase p(Prof::DirectLighting);
    // Randomly choose a single light to sample, _light_
    int nLights = scene.lights_size;
    if (nLights == 0) return Spectrum(0.f);
    int lightNum;
    Float lightPdf;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        lightNum = min((int)(sampler.Get1D() * nLights), nLights - 1);
        lightPdf = Float(1) / nLights;
    }
    Light const *light = scene.lights[lightNum];
    Point2f uLight = sampler.Get2D();
    Point2f uScattering = sampler.Get2D();
    return EstimateDirect(it, uScattering, *light, uLight,
                          scene, sampler, arena, handleMedia) / lightPdf;
}

__both__
bool BVHAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    // printf("nodes[%p]\n", nodes);
    // if (!nodes) return false;

    // // ProfilePhase p(Prof::AccelIntersect);
    // bool hit = false;
    // Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    // int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // // Follow ray through BVH nodes to find primitive intersections
    // int toVisitOffset = 0, currentNodeIndex = 0;
    // int nodesToVisit[64];
    // while (true) {
    //     const LinearBVHNode *node = &nodes[currentNodeIndex];
    //     // Check ray against BVH node
    //     if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
    //         if (node->nPrimitives > 0) {
    //             // Intersect ray with primitives in leaf BVH node
    //             for (int i = 0; i < node->nPrimitives; ++i)
    //                 if (primitives[node->primitivesOffset + i]->Intersect(
    //                         ray, isect))
    //                     hit = true;
    //             if (toVisitOffset == 0) break;
    //             currentNodeIndex = nodesToVisit[--toVisitOffset];
    //         } else {
    //             // Put far BVH node on _nodesToVisit_ stack, advance to near
    //             // node
    //             if (dirIsNeg[node->axis]) {
    //                 nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
    //                 currentNodeIndex = node->secondChildOffset;
    //             } else {
    //                 nodesToVisit[toVisitOffset++] = node->secondChildOffset;
    //                 currentNodeIndex = currentNodeIndex + 1;
    //             }
    //         }
    //     } else {
    //         if (toVisitOffset == 0) break;
    //         currentNodeIndex = nodesToVisit[--toVisitOffset];
    //     }
    // }
    //return hit;
    return false;
}


// scene
__both__
bool Scene::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    // ++nIntersectionTests;
    assert(ray.d != Vector3f(0,0,0));
    printf("ray.d[%p]\n", &(ray.d));
    printf("aggregate[%p]\n", aggregate);
    //return false;
    //printf("aggregate->Intersect[%p]\n", &(aggregate->Intersect));
    //return aggregate->Intersect(ray, isect);
    return false;
}

// interaction
__both__
Spectrum SurfaceInteraction::Le(const Vector3f &w) const {
    const AreaLight *area = primitive->GetAreaLight();
    return area ? area->L(*this, w) : Spectrum(0.f);
}





__global__
void LiKernel(Spectrum* Ls, PathIntegrator* integrator,
              const RayDifferential* rays, const Float* rayWeights,
              const Scene &scene, Sampler** tileSamplers,
              MemoryArena &arena) {
    int sampleNum = threadIdx.x;
    Ls[sampleNum] = 0.f;
    
    if (rayWeights[sampleNum] > 0){
        printf("integrator[%p], scene[%p], tileSampler[%p], arena[%p]\n", 
              integrator, &scene, tileSamplers[sampleNum], &arena);
        // printf("integrator.Li[%p]\n", &integrator->Li(rays[sampleNum], scene, *tileSamplers[sampleNum], arena, 0));
        Ls[sampleNum] = integrator->Li(rays[sampleNum], scene, *tileSamplers[sampleNum], arena, 0);
    }
}

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

//__noinline__
void Render(Integrator *i, Scene *s){
    PathIntegrator* integrator = dynamic_cast<PathIntegrator*>(i);
    Scene &scene = *s;
    printf("hereherehere\n");
    //CallLiKernel();

    integrator->Preprocess(scene, *(integrator->sampler));
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    Bounds2i sampleBounds = integrator->camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);
    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            void* arenaPtr;
            cudaMallocManaged(&arenaPtr, sizeof(MemoryArena));
            cudaDeviceSynchronize();
            MemoryArena* arena = new(arenaPtr) MemoryArena;

            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                integrator->camera->film->GetFilmTile(tileBounds);

            // Creating unified objects for calculations
            int samplesPerPixel = integrator->sampler->samplesPerPixel;
            void* tileSamplersPtr;
            Sampler** tileSamplers;
            CameraSample* cameraSamples;
            RayDifferential* rays;
            Float* rayWeights;
            Spectrum* Ls;
            cudaMallocManaged(&tileSamplersPtr, sizeof(*tileSamplers)
                              * samplesPerPixel);
            cudaMallocManaged(&cameraSamples, sizeof(CameraSample) 
                              * samplesPerPixel);
            cudaMallocManaged(&rays,          sizeof(RayDifferential) 
                              * samplesPerPixel); 
            cudaMallocManaged(&rayWeights,    sizeof(Float) 
                              * samplesPerPixel); 
            cudaMallocManaged(&Ls,            sizeof(Spectrum) 
                              * samplesPerPixel); 
            
            cudaDeviceSynchronize();
            new(cameraSamples) CameraSample[samplesPerPixel];
            new(rays)          RayDifferential[samplesPerPixel];
            new(rayWeights)    Float[samplesPerPixel];
            new(Ls)            Spectrum[samplesPerPixel];

            tileSamplers = (Sampler**)tileSamplersPtr;

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                if (!InsideExclusive(pixel, integrator->pixelBounds))
                    continue;

                for(int64_t sampleNum = 0; 
                    sampleNum < samplesPerPixel; 
                    sampleNum++) {

                    tileSamplers[sampleNum] = integrator->sampler->Clone(seed);
                    {
                        ProfilePhase pp(Prof::StartPixel);
                        tileSamplers[sampleNum]->StartPixel(pixel);
                    }
                    tileSamplers[sampleNum]->SetSampleNumber(sampleNum);
                
                    // Initialize _CameraSample_ for current sample
                    cameraSamples[sampleNum] =
                        tileSamplers[sampleNum]->GetCameraSample(pixel);

                    // Generate camera ray for current sample
                    
                    rayWeights[sampleNum] =
                        integrator->camera->GenerateRayDifferential(cameraSamples[sampleNum], &rays[sampleNum]);
                    rays[sampleNum].ScaleDifferentials(
                        1 / std::sqrt((Float)tileSamplers[sampleNum]->samplesPerPixel));
                    ++nCameraRays;
                }
                
                // Evaluate radiance along camera ray
                // CallLiKernel(s);
                
                cudaDeviceSynchronize();
                LiKernel<<<1, samplesPerPixel>>>
                    (Ls, integrator, rays, rayWeights, scene, tileSamplers, *arena);
                cudaDeviceSynchronize();

                // LOG(ERROR) << cudaGetErrorString(cudaGetLastError()) << std::endl;

                for(int64_t sampleNum = 0; 
                    sampleNum < samplesPerPixel; 
                    sampleNum++) {
                    
                    // Ls[sampleNum] = 0.f;
                    // if (rayWeights[sampleNum] > 0) 
                    //     Ls[sampleNum] = 
                    //         Li(rays[sampleNum], scene, *tileSamplers[sampleNum], arena);

                    // Issue warning if unexpected radiance value returned
                    // if (Ls[sampleNum] != Spectrum(20.0)) {
                    //     // LOG(ERROR) << "kernel may not be running\n";
                    // }

                    if (Ls[sampleNum].HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSamplers[sampleNum]->CurrentSampleNumber());
                        Ls[sampleNum] = Spectrum(0.f);
                    } else if (Ls[sampleNum].y() < -1e-5) {
                        LOG(ERROR) << StringPrintf(
                            "Negative luminance value, %f, returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            Ls[sampleNum].y(), pixel.x, pixel.y,
                            (int)tileSamplers[sampleNum]->CurrentSampleNumber());
                        Ls[sampleNum] = Spectrum(0.f);
                    } else if (std::isinf(Ls[sampleNum].y())) {
                          LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSamplers[sampleNum]->CurrentSampleNumber());
                        Ls[sampleNum] = Spectrum(0.f);
                    }
                    VLOG(1) << "Camera sample: " << cameraSamples[sampleNum] 
                            << " -> ray: " << rays[sampleNum]
                            << " -> Ls[sampleNum] = " << Ls[sampleNum];

                    // Add camera ray's contribution to image
                    filmTile->AddSample(cameraSamples[sampleNum].pFilm, 
                                        Ls[sampleNum], rayWeights[sampleNum]);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena->Reset();
                }
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            integrator->camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
            cudaFree(cameraSamples);
            cudaFree(rays);
            cudaFree(rayWeights);
            cudaFree(Ls);
        }, nTiles);
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
    integrator->camera->film->WriteImage();
}

Spectrum BSDF::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
                        const Point2f &u, Float *pdf, BxDFType type,
                        BxDFType *sampledType) const {
    // ProfilePhase pp(Prof::BSDFSampling);
    // Choose which _BxDF_ to sample
    int matchingComps = NumComponents(type);
    if (matchingComps == 0) {
        *pdf = 0;
        if (sampledType) *sampledType = BxDFType(0);
        return Spectrum(0);
    }
    int comp =
        min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

    // Get _BxDF_ pointer for chosen component
    BxDF *bxdf = nullptr;
    int count = comp;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
            bxdf = bxdfs[i];
            break;
        }
    // CHECK(bxdf != nullptr);
    // VLOG(2) << "BSDF::Sample_f chose comp = " << comp << " / matching = " <<
        // matchingComps << ", bxdf: " << bxdf->ToString();

    // Remap _BxDF_ sample _u_ to $[0,1)^2$
    Point2f uRemapped(min(u[0] * matchingComps - comp, OneMinusEpsilon),
                      u[1]);

    // Sample chosen _BxDF_
    Vector3f wi, wo = WorldToLocal(woWorld);
    if (wo.z == 0) return 0.;
    *pdf = 0;
    if (sampledType) *sampledType = bxdf->type;
    Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
    // VLOG(2) << "For wo = " << wo << ", sampled f = " << f << ", pdf = "
    //         << *pdf << ", ratio = " << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.))
    //         << ", wi = " << wi;
    if (*pdf == 0) {
        if (sampledType) *sampledType = BxDFType(0);
        return 0;
    }
    *wiWorld = LocalToWorld(wi);

    // Compute overall PDF with all matching _BxDF_s
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type))
                *pdf += bxdfs[i]->Pdf(wo, wi);
    if (matchingComps > 1) *pdf /= matchingComps;

    // Compute value of BSDF for sampled direction
    if (!(bxdf->type & BSDF_SPECULAR)) {
        bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
        f = 0.;
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(type) &&
                ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                 (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
                f += bxdfs[i]->f(wo, wi);
    }
    // VLOG(2) << "Overall f = " << f << ", pdf = " << *pdf << ", ratio = "
    //         << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.));
    return f;
}

Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW,
                 BxDFType flags) const {
    // ProfilePhase pp(Prof::BSDFEvaluation);
    Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
    if (wo.z == 0) return 0.;
    bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
    Spectrum f(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags) &&
            ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
             (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
            f += bxdfs[i]->f(wo, wi);
    return f;
}

Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
                BxDFType flags) const {
    // ProfilePhase pp(Prof::BSDFPdf);
    if (nBxDFs == 0.f) return 0.f;
    Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
    if (wo.z == 0) return 0.;
    Float pdf = 0.f;
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
    Float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
    return v;
}

__both__
bool VisibilityTester::Unoccluded(const Scene &scene) const {
    return !scene.IntersectP(p0.SpawnRayTo(p1));
}

__both__
void SurfaceInteraction::ComputeScatteringFunctions(const RayDifferential &ray,
                                                    MemoryArena &arena,
                                                    bool allowMultipleLobes,
                                                    TransportMode mode) {
    ComputeDifferentials(ray);
    primitive->ComputeScatteringFunctions(this, arena, mode,
                                          allowMultipleLobes);
}

__both__
Spectrum Light::Le(const RayDifferential &ray) const { return Spectrum(0.f); }

Spectrum BxDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                        Float *pdf, BxDFType *sampledType) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere(u);
    if (wo.z < 0) wi->z *= -1;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
    Spectrum r(0.);
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hd}$
        Vector3f wi;
        Float pdf = 0;
        Spectrum f = Sample_f(w, &wi, u[i], &pdf);
        if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
    }
    return r / nSamples;
}

Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
    Spectrum r(0.f);
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hh}$
        Vector3f wo, wi;
        wo = UniformSampleHemisphere(u1[i]);
        Float pdfo = UniformHemispherePdf(), pdfi = 0;
        Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
        if (pdfi > 0)
            r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
    }
    return r / (Pi * nSamples);
}

// BxDF Method Definitions
Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    return R * InvPi;
}

__both__
void Material::Bump(Texture<Float> const * d,
                    SurfaceInteraction *si) {
    // Compute offset positions and evaluate displacement texture
    SurfaceInteraction siEval = *si;

    // Shift _siEval_ _du_ in the $u$ direction
    Float du = .5f * (pbrt::math::abs(si->dudx) + pbrt::math::abs(si->dudy));
    // The most common reason for du to be zero is for ray that start from
    // light sources, where no differentials are available. In this case,
    // we try to choose a small enough du so that we still get a decently
    // accurate bump value.
    if (du == 0) du = .0005f;
    siEval.p = si->p + du * si->shading.dpdu;
    siEval.uv = si->uv + Vector2f(du, 0.f);
    siEval.n = Normalize((Normal3f)Cross(si->shading.dpdu, si->shading.dpdv) +
                         du * si->dndu);
    Float uDisplace = d->Evaluate(siEval);

    // Shift _siEval_ _dv_ in the $v$ direction
    Float dv = .5f * (pbrt::math::abs(si->dvdx) + pbrt::math::abs(si->dvdy));
    if (dv == 0) dv = .0005f;
    siEval.p = si->p + dv * si->shading.dpdv;
    siEval.uv = si->uv + Vector2f(0.f, dv);
    siEval.n = Normalize((Normal3f)Cross(si->shading.dpdu, si->shading.dpdv) +
                         dv * si->dndv);
    Float vDisplace = d->Evaluate(siEval);
    Float displace = d->Evaluate(*si);

    // Compute bump-mapped differential geometry
    Vector3f dpdu = si->shading.dpdu +
                    (uDisplace - displace) / du * Vector3f(si->shading.n) +
                    displace * Vector3f(si->shading.dndu);
    Vector3f dpdv = si->shading.dpdv +
                    (vDisplace - displace) / dv * Vector3f(si->shading.n) +
                    displace * Vector3f(si->shading.dndv);
    si->SetShadingGeometry(dpdu, dpdv, si->shading.dndu, si->shading.dndv,
                           false);
}

Float GlobalSampler::Get1D() {
    // ProfilePhase _(Prof::GetSample);
    if (dimension >= arrayStartDim && dimension < arrayEndDim)
        dimension = arrayEndDim;
    return SampleDimension(intervalSampleIndex, dimension++);
}

Point2f GlobalSampler::Get2D() {
    // ProfilePhase _(Prof::GetSample);
    if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim)
        dimension = arrayEndDim;
    Point2f p(SampleDimension(intervalSampleIndex, dimension),
              SampleDimension(intervalSampleIndex, dimension + 1));
    dimension += 2;
    return p;
}

template <int base>
__both__
PBRT_NOINLINE static Float RadicalInverseSpecialized(uint64_t a) {
    const Float invBase = (Float)1 / (Float)base;
    uint64_t reversedDigits = 0;
    Float invBaseN = 1;
    while (a) {
        uint64_t next = a / base;
        uint64_t digit = a - next * base;
        reversedDigits = reversedDigits * base + digit;
        invBaseN *= invBase;
        a = next;
    }
    assert(reversedDigits * invBaseN < 1.00001);
    return min(reversedDigits * invBaseN, OneMinusEpsilon);
}

template <int base>
__both__
PBRT_NOINLINE static Float
ScrambledRadicalInverseSpecialized(const uint16_t *perm, uint64_t a) {
    const Float invBase = (Float)1 / (Float)base;
    uint64_t reversedDigits = 0;
    Float invBaseN = 1;
    while (a) {
        uint64_t next = a / base;
        uint64_t digit = a - next * base;
        // CHECK_LT(perm[digit], base);
        reversedDigits = reversedDigits * base + perm[digit];
        invBaseN *= invBase;
        a = next;
    }
    assert(invBaseN * (reversedDigits + invBase * perm[0] / (1 - invBase))
              < 1.00001);
    return min(
        invBaseN * (reversedDigits + invBase * perm[0] / (1 - invBase)),
        OneMinusEpsilon);
}

__both__
Float RadicalInverse(int baseIndex, uint64_t a) {
    switch (baseIndex) {
    case 0:
    // Compute base-2 radical inverse
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        return ReverseBits64(a) * 5.4210108624275222e-20;
#else
        return ReverseBits64(a) * 0x1p-64;
#endif
    case 1:
        return RadicalInverseSpecialized<3>(a);
    case 2:
        return RadicalInverseSpecialized<5>(a);
    case 3:
        return RadicalInverseSpecialized<7>(a);
    // Remainder of cases for _RadicalInverse()_
    case 4:
        return RadicalInverseSpecialized<11>(a);
    case 5:
        return RadicalInverseSpecialized<13>(a);
    case 6:
        return RadicalInverseSpecialized<17>(a);
    case 7:
        return RadicalInverseSpecialized<19>(a);
    case 8:
        return RadicalInverseSpecialized<23>(a);
    case 9:
        return RadicalInverseSpecialized<29>(a);
    case 10:
        return RadicalInverseSpecialized<31>(a);
    case 11:
        return RadicalInverseSpecialized<37>(a);
    case 12:
        return RadicalInverseSpecialized<41>(a);
    case 13:
        return RadicalInverseSpecialized<43>(a);
    case 14:
        return RadicalInverseSpecialized<47>(a);
    case 15:
        return RadicalInverseSpecialized<53>(a);
    case 16:
        return RadicalInverseSpecialized<59>(a);
    case 17:
        return RadicalInverseSpecialized<61>(a);
    case 18:
        return RadicalInverseSpecialized<67>(a);
    case 19:
        return RadicalInverseSpecialized<71>(a);
    case 20:
        return RadicalInverseSpecialized<73>(a);
    case 21:
        return RadicalInverseSpecialized<79>(a);
    case 22:
        return RadicalInverseSpecialized<83>(a);
    case 23:
        return RadicalInverseSpecialized<89>(a);
    case 24:
        return RadicalInverseSpecialized<97>(a);
    case 25:
        return RadicalInverseSpecialized<101>(a);
    case 26:
        return RadicalInverseSpecialized<103>(a);
    case 27:
        return RadicalInverseSpecialized<107>(a);
    case 28:
        return RadicalInverseSpecialized<109>(a);
    case 29:
        return RadicalInverseSpecialized<113>(a);
    case 30:
        return RadicalInverseSpecialized<127>(a);
    case 31:
        return RadicalInverseSpecialized<131>(a);
    case 32:
        return RadicalInverseSpecialized<137>(a);
    case 33:
        return RadicalInverseSpecialized<139>(a);
    case 34:
        return RadicalInverseSpecialized<149>(a);
    case 35:
        return RadicalInverseSpecialized<151>(a);
    case 36:
        return RadicalInverseSpecialized<157>(a);
    case 37:
        return RadicalInverseSpecialized<163>(a);
    case 38:
        return RadicalInverseSpecialized<167>(a);
    case 39:
        return RadicalInverseSpecialized<173>(a);
    case 40:
        return RadicalInverseSpecialized<179>(a);
    case 41:
        return RadicalInverseSpecialized<181>(a);
    case 42:
        return RadicalInverseSpecialized<191>(a);
    case 43:
        return RadicalInverseSpecialized<193>(a);
    case 44:
        return RadicalInverseSpecialized<197>(a);
    case 45:
        return RadicalInverseSpecialized<199>(a);
    case 46:
        return RadicalInverseSpecialized<211>(a);
    case 47:
        return RadicalInverseSpecialized<223>(a);
    case 48:
        return RadicalInverseSpecialized<227>(a);
    case 49:
        return RadicalInverseSpecialized<229>(a);
    case 50:
        return RadicalInverseSpecialized<233>(a);
    case 51:
        return RadicalInverseSpecialized<239>(a);
    case 52:
        return RadicalInverseSpecialized<241>(a);
    case 53:
        return RadicalInverseSpecialized<251>(a);
    case 54:
        return RadicalInverseSpecialized<257>(a);
    case 55:
        return RadicalInverseSpecialized<263>(a);
    case 56:
        return RadicalInverseSpecialized<269>(a);
    case 57:
        return RadicalInverseSpecialized<271>(a);
    case 58:
        return RadicalInverseSpecialized<277>(a);
    case 59:
        return RadicalInverseSpecialized<281>(a);
    case 60:
        return RadicalInverseSpecialized<283>(a);
    case 61:
        return RadicalInverseSpecialized<293>(a);
    case 62:
        return RadicalInverseSpecialized<307>(a);
    case 63:
        return RadicalInverseSpecialized<311>(a);
    case 64:
        return RadicalInverseSpecialized<313>(a);
    case 65:
        return RadicalInverseSpecialized<317>(a);
    case 66:
        return RadicalInverseSpecialized<331>(a);
    case 67:
        return RadicalInverseSpecialized<337>(a);
    case 68:
        return RadicalInverseSpecialized<347>(a);
    case 69:
        return RadicalInverseSpecialized<349>(a);
    case 70:
        return RadicalInverseSpecialized<353>(a);
    case 71:
        return RadicalInverseSpecialized<359>(a);
    case 72:
        return RadicalInverseSpecialized<367>(a);
    case 73:
        return RadicalInverseSpecialized<373>(a);
    case 74:
        return RadicalInverseSpecialized<379>(a);
    case 75:
        return RadicalInverseSpecialized<383>(a);
    case 76:
        return RadicalInverseSpecialized<389>(a);
    case 77:
        return RadicalInverseSpecialized<397>(a);
    case 78:
        return RadicalInverseSpecialized<401>(a);
    case 79:
        return RadicalInverseSpecialized<409>(a);
    case 80:
        return RadicalInverseSpecialized<419>(a);
    case 81:
        return RadicalInverseSpecialized<421>(a);
    case 82:
        return RadicalInverseSpecialized<431>(a);
    case 83:
        return RadicalInverseSpecialized<433>(a);
    case 84:
        return RadicalInverseSpecialized<439>(a);
    case 85:
        return RadicalInverseSpecialized<443>(a);
    case 86:
        return RadicalInverseSpecialized<449>(a);
    case 87:
        return RadicalInverseSpecialized<457>(a);
    case 88:
        return RadicalInverseSpecialized<461>(a);
    case 89:
        return RadicalInverseSpecialized<463>(a);
    case 90:
        return RadicalInverseSpecialized<467>(a);
    case 91:
        return RadicalInverseSpecialized<479>(a);
    case 92:
        return RadicalInverseSpecialized<487>(a);
    case 93:
        return RadicalInverseSpecialized<491>(a);
    case 94:
        return RadicalInverseSpecialized<499>(a);
    case 95:
        return RadicalInverseSpecialized<503>(a);
    case 96:
        return RadicalInverseSpecialized<509>(a);
    case 97:
        return RadicalInverseSpecialized<521>(a);
    case 98:
        return RadicalInverseSpecialized<523>(a);
    case 99:
        return RadicalInverseSpecialized<541>(a);
    case 100:
        return RadicalInverseSpecialized<547>(a);
    case 101:
        return RadicalInverseSpecialized<557>(a);
    case 102:
        return RadicalInverseSpecialized<563>(a);
    case 103:
        return RadicalInverseSpecialized<569>(a);
    case 104:
        return RadicalInverseSpecialized<571>(a);
    case 105:
        return RadicalInverseSpecialized<577>(a);
    case 106:
        return RadicalInverseSpecialized<587>(a);
    case 107:
        return RadicalInverseSpecialized<593>(a);
    case 108:
        return RadicalInverseSpecialized<599>(a);
    case 109:
        return RadicalInverseSpecialized<601>(a);
    case 110:
        return RadicalInverseSpecialized<607>(a);
    case 111:
        return RadicalInverseSpecialized<613>(a);
    case 112:
        return RadicalInverseSpecialized<617>(a);
    case 113:
        return RadicalInverseSpecialized<619>(a);
    case 114:
        return RadicalInverseSpecialized<631>(a);
    case 115:
        return RadicalInverseSpecialized<641>(a);
    case 116:
        return RadicalInverseSpecialized<643>(a);
    case 117:
        return RadicalInverseSpecialized<647>(a);
    case 118:
        return RadicalInverseSpecialized<653>(a);
    case 119:
        return RadicalInverseSpecialized<659>(a);
    case 120:
        return RadicalInverseSpecialized<661>(a);
    case 121:
        return RadicalInverseSpecialized<673>(a);
    case 122:
        return RadicalInverseSpecialized<677>(a);
    case 123:
        return RadicalInverseSpecialized<683>(a);
    case 124:
        return RadicalInverseSpecialized<691>(a);
    case 125:
        return RadicalInverseSpecialized<701>(a);
    case 126:
        return RadicalInverseSpecialized<709>(a);
    case 127:
        return RadicalInverseSpecialized<719>(a);
    case 128:
        return RadicalInverseSpecialized<727>(a);
    case 129:
        return RadicalInverseSpecialized<733>(a);
    case 130:
        return RadicalInverseSpecialized<739>(a);
    case 131:
        return RadicalInverseSpecialized<743>(a);
    case 132:
        return RadicalInverseSpecialized<751>(a);
    case 133:
        return RadicalInverseSpecialized<757>(a);
    case 134:
        return RadicalInverseSpecialized<761>(a);
    case 135:
        return RadicalInverseSpecialized<769>(a);
    case 136:
        return RadicalInverseSpecialized<773>(a);
    case 137:
        return RadicalInverseSpecialized<787>(a);
    case 138:
        return RadicalInverseSpecialized<797>(a);
    case 139:
        return RadicalInverseSpecialized<809>(a);
    case 140:
        return RadicalInverseSpecialized<811>(a);
    case 141:
        return RadicalInverseSpecialized<821>(a);
    case 142:
        return RadicalInverseSpecialized<823>(a);
    case 143:
        return RadicalInverseSpecialized<827>(a);
    case 144:
        return RadicalInverseSpecialized<829>(a);
    case 145:
        return RadicalInverseSpecialized<839>(a);
    case 146:
        return RadicalInverseSpecialized<853>(a);
    case 147:
        return RadicalInverseSpecialized<857>(a);
    case 148:
        return RadicalInverseSpecialized<859>(a);
    case 149:
        return RadicalInverseSpecialized<863>(a);
    case 150:
        return RadicalInverseSpecialized<877>(a);
    case 151:
        return RadicalInverseSpecialized<881>(a);
    case 152:
        return RadicalInverseSpecialized<883>(a);
    case 153:
        return RadicalInverseSpecialized<887>(a);
    case 154:
        return RadicalInverseSpecialized<907>(a);
    case 155:
        return RadicalInverseSpecialized<911>(a);
    case 156:
        return RadicalInverseSpecialized<919>(a);
    case 157:
        return RadicalInverseSpecialized<929>(a);
    case 158:
        return RadicalInverseSpecialized<937>(a);
    case 159:
        return RadicalInverseSpecialized<941>(a);
    case 160:
        return RadicalInverseSpecialized<947>(a);
    case 161:
        return RadicalInverseSpecialized<953>(a);
    case 162:
        return RadicalInverseSpecialized<967>(a);
    case 163:
        return RadicalInverseSpecialized<971>(a);
    case 164:
        return RadicalInverseSpecialized<977>(a);
    case 165:
        return RadicalInverseSpecialized<983>(a);
    case 166:
        return RadicalInverseSpecialized<991>(a);
    case 167:
        return RadicalInverseSpecialized<997>(a);
    case 168:
        return RadicalInverseSpecialized<1009>(a);
    case 169:
        return RadicalInverseSpecialized<1013>(a);
    case 170:
        return RadicalInverseSpecialized<1019>(a);
    case 171:
        return RadicalInverseSpecialized<1021>(a);
    case 172:
        return RadicalInverseSpecialized<1031>(a);
    case 173:
        return RadicalInverseSpecialized<1033>(a);
    case 174:
        return RadicalInverseSpecialized<1039>(a);
    case 175:
        return RadicalInverseSpecialized<1049>(a);
    case 176:
        return RadicalInverseSpecialized<1051>(a);
    case 177:
        return RadicalInverseSpecialized<1061>(a);
    case 178:
        return RadicalInverseSpecialized<1063>(a);
    case 179:
        return RadicalInverseSpecialized<1069>(a);
    case 180:
        return RadicalInverseSpecialized<1087>(a);
    case 181:
        return RadicalInverseSpecialized<1091>(a);
    case 182:
        return RadicalInverseSpecialized<1093>(a);
    case 183:
        return RadicalInverseSpecialized<1097>(a);
    case 184:
        return RadicalInverseSpecialized<1103>(a);
    case 185:
        return RadicalInverseSpecialized<1109>(a);
    case 186:
        return RadicalInverseSpecialized<1117>(a);
    case 187:
        return RadicalInverseSpecialized<1123>(a);
    case 188:
        return RadicalInverseSpecialized<1129>(a);
    case 189:
        return RadicalInverseSpecialized<1151>(a);
    case 190:
        return RadicalInverseSpecialized<1153>(a);
    case 191:
        return RadicalInverseSpecialized<1163>(a);
    case 192:
        return RadicalInverseSpecialized<1171>(a);
    case 193:
        return RadicalInverseSpecialized<1181>(a);
    case 194:
        return RadicalInverseSpecialized<1187>(a);
    case 195:
        return RadicalInverseSpecialized<1193>(a);
    case 196:
        return RadicalInverseSpecialized<1201>(a);
    case 197:
        return RadicalInverseSpecialized<1213>(a);
    case 198:
        return RadicalInverseSpecialized<1217>(a);
    case 199:
        return RadicalInverseSpecialized<1223>(a);
    case 200:
        return RadicalInverseSpecialized<1229>(a);
    case 201:
        return RadicalInverseSpecialized<1231>(a);
    case 202:
        return RadicalInverseSpecialized<1237>(a);
    case 203:
        return RadicalInverseSpecialized<1249>(a);
    case 204:
        return RadicalInverseSpecialized<1259>(a);
    case 205:
        return RadicalInverseSpecialized<1277>(a);
    case 206:
        return RadicalInverseSpecialized<1279>(a);
    case 207:
        return RadicalInverseSpecialized<1283>(a);
    case 208:
        return RadicalInverseSpecialized<1289>(a);
    case 209:
        return RadicalInverseSpecialized<1291>(a);
    case 210:
        return RadicalInverseSpecialized<1297>(a);
    case 211:
        return RadicalInverseSpecialized<1301>(a);
    case 212:
        return RadicalInverseSpecialized<1303>(a);
    case 213:
        return RadicalInverseSpecialized<1307>(a);
    case 214:
        return RadicalInverseSpecialized<1319>(a);
    case 215:
        return RadicalInverseSpecialized<1321>(a);
    case 216:
        return RadicalInverseSpecialized<1327>(a);
    case 217:
        return RadicalInverseSpecialized<1361>(a);
    case 218:
        return RadicalInverseSpecialized<1367>(a);
    case 219:
        return RadicalInverseSpecialized<1373>(a);
    case 220:
        return RadicalInverseSpecialized<1381>(a);
    case 221:
        return RadicalInverseSpecialized<1399>(a);
    case 222:
        return RadicalInverseSpecialized<1409>(a);
    case 223:
        return RadicalInverseSpecialized<1423>(a);
    case 224:
        return RadicalInverseSpecialized<1427>(a);
    case 225:
        return RadicalInverseSpecialized<1429>(a);
    case 226:
        return RadicalInverseSpecialized<1433>(a);
    case 227:
        return RadicalInverseSpecialized<1439>(a);
    case 228:
        return RadicalInverseSpecialized<1447>(a);
    case 229:
        return RadicalInverseSpecialized<1451>(a);
    case 230:
        return RadicalInverseSpecialized<1453>(a);
    case 231:
        return RadicalInverseSpecialized<1459>(a);
    case 232:
        return RadicalInverseSpecialized<1471>(a);
    case 233:
        return RadicalInverseSpecialized<1481>(a);
    case 234:
        return RadicalInverseSpecialized<1483>(a);
    case 235:
        return RadicalInverseSpecialized<1487>(a);
    case 236:
        return RadicalInverseSpecialized<1489>(a);
    case 237:
        return RadicalInverseSpecialized<1493>(a);
    case 238:
        return RadicalInverseSpecialized<1499>(a);
    case 239:
        return RadicalInverseSpecialized<1511>(a);
    case 240:
        return RadicalInverseSpecialized<1523>(a);
    case 241:
        return RadicalInverseSpecialized<1531>(a);
    case 242:
        return RadicalInverseSpecialized<1543>(a);
    case 243:
        return RadicalInverseSpecialized<1549>(a);
    case 244:
        return RadicalInverseSpecialized<1553>(a);
    case 245:
        return RadicalInverseSpecialized<1559>(a);
    case 246:
        return RadicalInverseSpecialized<1567>(a);
    case 247:
        return RadicalInverseSpecialized<1571>(a);
    case 248:
        return RadicalInverseSpecialized<1579>(a);
    case 249:
        return RadicalInverseSpecialized<1583>(a);
    case 250:
        return RadicalInverseSpecialized<1597>(a);
    case 251:
        return RadicalInverseSpecialized<1601>(a);
    case 252:
        return RadicalInverseSpecialized<1607>(a);
    case 253:
        return RadicalInverseSpecialized<1609>(a);
    case 254:
        return RadicalInverseSpecialized<1613>(a);
    case 255:
        return RadicalInverseSpecialized<1619>(a);
    case 256:
        return RadicalInverseSpecialized<1621>(a);
    case 257:
        return RadicalInverseSpecialized<1627>(a);
    case 258:
        return RadicalInverseSpecialized<1637>(a);
    case 259:
        return RadicalInverseSpecialized<1657>(a);
    case 260:
        return RadicalInverseSpecialized<1663>(a);
    case 261:
        return RadicalInverseSpecialized<1667>(a);
    case 262:
        return RadicalInverseSpecialized<1669>(a);
    case 263:
        return RadicalInverseSpecialized<1693>(a);
    case 264:
        return RadicalInverseSpecialized<1697>(a);
    case 265:
        return RadicalInverseSpecialized<1699>(a);
    case 266:
        return RadicalInverseSpecialized<1709>(a);
    case 267:
        return RadicalInverseSpecialized<1721>(a);
    case 268:
        return RadicalInverseSpecialized<1723>(a);
    case 269:
        return RadicalInverseSpecialized<1733>(a);
    case 270:
        return RadicalInverseSpecialized<1741>(a);
    case 271:
        return RadicalInverseSpecialized<1747>(a);
    case 272:
        return RadicalInverseSpecialized<1753>(a);
    case 273:
        return RadicalInverseSpecialized<1759>(a);
    case 274:
        return RadicalInverseSpecialized<1777>(a);
    case 275:
        return RadicalInverseSpecialized<1783>(a);
    case 276:
        return RadicalInverseSpecialized<1787>(a);
    case 277:
        return RadicalInverseSpecialized<1789>(a);
    case 278:
        return RadicalInverseSpecialized<1801>(a);
    case 279:
        return RadicalInverseSpecialized<1811>(a);
    case 280:
        return RadicalInverseSpecialized<1823>(a);
    case 281:
        return RadicalInverseSpecialized<1831>(a);
    case 282:
        return RadicalInverseSpecialized<1847>(a);
    case 283:
        return RadicalInverseSpecialized<1861>(a);
    case 284:
        return RadicalInverseSpecialized<1867>(a);
    case 285:
        return RadicalInverseSpecialized<1871>(a);
    case 286:
        return RadicalInverseSpecialized<1873>(a);
    case 287:
        return RadicalInverseSpecialized<1877>(a);
    case 288:
        return RadicalInverseSpecialized<1879>(a);
    case 289:
        return RadicalInverseSpecialized<1889>(a);
    case 290:
        return RadicalInverseSpecialized<1901>(a);
    case 291:
        return RadicalInverseSpecialized<1907>(a);
    case 292:
        return RadicalInverseSpecialized<1913>(a);
    case 293:
        return RadicalInverseSpecialized<1931>(a);
    case 294:
        return RadicalInverseSpecialized<1933>(a);
    case 295:
        return RadicalInverseSpecialized<1949>(a);
    case 296:
        return RadicalInverseSpecialized<1951>(a);
    case 297:
        return RadicalInverseSpecialized<1973>(a);
    case 298:
        return RadicalInverseSpecialized<1979>(a);
    case 299:
        return RadicalInverseSpecialized<1987>(a);
    case 300:
        return RadicalInverseSpecialized<1993>(a);
    case 301:
        return RadicalInverseSpecialized<1997>(a);
    case 302:
        return RadicalInverseSpecialized<1999>(a);
    case 303:
        return RadicalInverseSpecialized<2003>(a);
    case 304:
        return RadicalInverseSpecialized<2011>(a);
    case 305:
        return RadicalInverseSpecialized<2017>(a);
    case 306:
        return RadicalInverseSpecialized<2027>(a);
    case 307:
        return RadicalInverseSpecialized<2029>(a);
    case 308:
        return RadicalInverseSpecialized<2039>(a);
    case 309:
        return RadicalInverseSpecialized<2053>(a);
    case 310:
        return RadicalInverseSpecialized<2063>(a);
    case 311:
        return RadicalInverseSpecialized<2069>(a);
    case 312:
        return RadicalInverseSpecialized<2081>(a);
    case 313:
        return RadicalInverseSpecialized<2083>(a);
    case 314:
        return RadicalInverseSpecialized<2087>(a);
    case 315:
        return RadicalInverseSpecialized<2089>(a);
    case 316:
        return RadicalInverseSpecialized<2099>(a);
    case 317:
        return RadicalInverseSpecialized<2111>(a);
    case 318:
        return RadicalInverseSpecialized<2113>(a);
    case 319:
        return RadicalInverseSpecialized<2129>(a);
    case 320:
        return RadicalInverseSpecialized<2131>(a);
    case 321:
        return RadicalInverseSpecialized<2137>(a);
    case 322:
        return RadicalInverseSpecialized<2141>(a);
    case 323:
        return RadicalInverseSpecialized<2143>(a);
    case 324:
        return RadicalInverseSpecialized<2153>(a);
    case 325:
        return RadicalInverseSpecialized<2161>(a);
    case 326:
        return RadicalInverseSpecialized<2179>(a);
    case 327:
        return RadicalInverseSpecialized<2203>(a);
    case 328:
        return RadicalInverseSpecialized<2207>(a);
    case 329:
        return RadicalInverseSpecialized<2213>(a);
    case 330:
        return RadicalInverseSpecialized<2221>(a);
    case 331:
        return RadicalInverseSpecialized<2237>(a);
    case 332:
        return RadicalInverseSpecialized<2239>(a);
    case 333:
        return RadicalInverseSpecialized<2243>(a);
    case 334:
        return RadicalInverseSpecialized<2251>(a);
    case 335:
        return RadicalInverseSpecialized<2267>(a);
    case 336:
        return RadicalInverseSpecialized<2269>(a);
    case 337:
        return RadicalInverseSpecialized<2273>(a);
    case 338:
        return RadicalInverseSpecialized<2281>(a);
    case 339:
        return RadicalInverseSpecialized<2287>(a);
    case 340:
        return RadicalInverseSpecialized<2293>(a);
    case 341:
        return RadicalInverseSpecialized<2297>(a);
    case 342:
        return RadicalInverseSpecialized<2309>(a);
    case 343:
        return RadicalInverseSpecialized<2311>(a);
    case 344:
        return RadicalInverseSpecialized<2333>(a);
    case 345:
        return RadicalInverseSpecialized<2339>(a);
    case 346:
        return RadicalInverseSpecialized<2341>(a);
    case 347:
        return RadicalInverseSpecialized<2347>(a);
    case 348:
        return RadicalInverseSpecialized<2351>(a);
    case 349:
        return RadicalInverseSpecialized<2357>(a);
    case 350:
        return RadicalInverseSpecialized<2371>(a);
    case 351:
        return RadicalInverseSpecialized<2377>(a);
    case 352:
        return RadicalInverseSpecialized<2381>(a);
    case 353:
        return RadicalInverseSpecialized<2383>(a);
    case 354:
        return RadicalInverseSpecialized<2389>(a);
    case 355:
        return RadicalInverseSpecialized<2393>(a);
    case 356:
        return RadicalInverseSpecialized<2399>(a);
    case 357:
        return RadicalInverseSpecialized<2411>(a);
    case 358:
        return RadicalInverseSpecialized<2417>(a);
    case 359:
        return RadicalInverseSpecialized<2423>(a);
    case 360:
        return RadicalInverseSpecialized<2437>(a);
    case 361:
        return RadicalInverseSpecialized<2441>(a);
    case 362:
        return RadicalInverseSpecialized<2447>(a);
    case 363:
        return RadicalInverseSpecialized<2459>(a);
    case 364:
        return RadicalInverseSpecialized<2467>(a);
    case 365:
        return RadicalInverseSpecialized<2473>(a);
    case 366:
        return RadicalInverseSpecialized<2477>(a);
    case 367:
        return RadicalInverseSpecialized<2503>(a);
    case 368:
        return RadicalInverseSpecialized<2521>(a);
    case 369:
        return RadicalInverseSpecialized<2531>(a);
    case 370:
        return RadicalInverseSpecialized<2539>(a);
    case 371:
        return RadicalInverseSpecialized<2543>(a);
    case 372:
        return RadicalInverseSpecialized<2549>(a);
    case 373:
        return RadicalInverseSpecialized<2551>(a);
    case 374:
        return RadicalInverseSpecialized<2557>(a);
    case 375:
        return RadicalInverseSpecialized<2579>(a);
    case 376:
        return RadicalInverseSpecialized<2591>(a);
    case 377:
        return RadicalInverseSpecialized<2593>(a);
    case 378:
        return RadicalInverseSpecialized<2609>(a);
    case 379:
        return RadicalInverseSpecialized<2617>(a);
    case 380:
        return RadicalInverseSpecialized<2621>(a);
    case 381:
        return RadicalInverseSpecialized<2633>(a);
    case 382:
        return RadicalInverseSpecialized<2647>(a);
    case 383:
        return RadicalInverseSpecialized<2657>(a);
    case 384:
        return RadicalInverseSpecialized<2659>(a);
    case 385:
        return RadicalInverseSpecialized<2663>(a);
    case 386:
        return RadicalInverseSpecialized<2671>(a);
    case 387:
        return RadicalInverseSpecialized<2677>(a);
    case 388:
        return RadicalInverseSpecialized<2683>(a);
    case 389:
        return RadicalInverseSpecialized<2687>(a);
    case 390:
        return RadicalInverseSpecialized<2689>(a);
    case 391:
        return RadicalInverseSpecialized<2693>(a);
    case 392:
        return RadicalInverseSpecialized<2699>(a);
    case 393:
        return RadicalInverseSpecialized<2707>(a);
    case 394:
        return RadicalInverseSpecialized<2711>(a);
    case 395:
        return RadicalInverseSpecialized<2713>(a);
    case 396:
        return RadicalInverseSpecialized<2719>(a);
    case 397:
        return RadicalInverseSpecialized<2729>(a);
    case 398:
        return RadicalInverseSpecialized<2731>(a);
    case 399:
        return RadicalInverseSpecialized<2741>(a);
    case 400:
        return RadicalInverseSpecialized<2749>(a);
    case 401:
        return RadicalInverseSpecialized<2753>(a);
    case 402:
        return RadicalInverseSpecialized<2767>(a);
    case 403:
        return RadicalInverseSpecialized<2777>(a);
    case 404:
        return RadicalInverseSpecialized<2789>(a);
    case 405:
        return RadicalInverseSpecialized<2791>(a);
    case 406:
        return RadicalInverseSpecialized<2797>(a);
    case 407:
        return RadicalInverseSpecialized<2801>(a);
    case 408:
        return RadicalInverseSpecialized<2803>(a);
    case 409:
        return RadicalInverseSpecialized<2819>(a);
    case 410:
        return RadicalInverseSpecialized<2833>(a);
    case 411:
        return RadicalInverseSpecialized<2837>(a);
    case 412:
        return RadicalInverseSpecialized<2843>(a);
    case 413:
        return RadicalInverseSpecialized<2851>(a);
    case 414:
        return RadicalInverseSpecialized<2857>(a);
    case 415:
        return RadicalInverseSpecialized<2861>(a);
    case 416:
        return RadicalInverseSpecialized<2879>(a);
    case 417:
        return RadicalInverseSpecialized<2887>(a);
    case 418:
        return RadicalInverseSpecialized<2897>(a);
    case 419:
        return RadicalInverseSpecialized<2903>(a);
    case 420:
        return RadicalInverseSpecialized<2909>(a);
    case 421:
        return RadicalInverseSpecialized<2917>(a);
    case 422:
        return RadicalInverseSpecialized<2927>(a);
    case 423:
        return RadicalInverseSpecialized<2939>(a);
    case 424:
        return RadicalInverseSpecialized<2953>(a);
    case 425:
        return RadicalInverseSpecialized<2957>(a);
    case 426:
        return RadicalInverseSpecialized<2963>(a);
    case 427:
        return RadicalInverseSpecialized<2969>(a);
    case 428:
        return RadicalInverseSpecialized<2971>(a);
    case 429:
        return RadicalInverseSpecialized<2999>(a);
    case 430:
        return RadicalInverseSpecialized<3001>(a);
    case 431:
        return RadicalInverseSpecialized<3011>(a);
    case 432:
        return RadicalInverseSpecialized<3019>(a);
    case 433:
        return RadicalInverseSpecialized<3023>(a);
    case 434:
        return RadicalInverseSpecialized<3037>(a);
    case 435:
        return RadicalInverseSpecialized<3041>(a);
    case 436:
        return RadicalInverseSpecialized<3049>(a);
    case 437:
        return RadicalInverseSpecialized<3061>(a);
    case 438:
        return RadicalInverseSpecialized<3067>(a);
    case 439:
        return RadicalInverseSpecialized<3079>(a);
    case 440:
        return RadicalInverseSpecialized<3083>(a);
    case 441:
        return RadicalInverseSpecialized<3089>(a);
    case 442:
        return RadicalInverseSpecialized<3109>(a);
    case 443:
        return RadicalInverseSpecialized<3119>(a);
    case 444:
        return RadicalInverseSpecialized<3121>(a);
    case 445:
        return RadicalInverseSpecialized<3137>(a);
    case 446:
        return RadicalInverseSpecialized<3163>(a);
    case 447:
        return RadicalInverseSpecialized<3167>(a);
    case 448:
        return RadicalInverseSpecialized<3169>(a);
    case 449:
        return RadicalInverseSpecialized<3181>(a);
    case 450:
        return RadicalInverseSpecialized<3187>(a);
    case 451:
        return RadicalInverseSpecialized<3191>(a);
    case 452:
        return RadicalInverseSpecialized<3203>(a);
    case 453:
        return RadicalInverseSpecialized<3209>(a);
    case 454:
        return RadicalInverseSpecialized<3217>(a);
    case 455:
        return RadicalInverseSpecialized<3221>(a);
    case 456:
        return RadicalInverseSpecialized<3229>(a);
    case 457:
        return RadicalInverseSpecialized<3251>(a);
    case 458:
        return RadicalInverseSpecialized<3253>(a);
    case 459:
        return RadicalInverseSpecialized<3257>(a);
    case 460:
        return RadicalInverseSpecialized<3259>(a);
    case 461:
        return RadicalInverseSpecialized<3271>(a);
    case 462:
        return RadicalInverseSpecialized<3299>(a);
    case 463:
        return RadicalInverseSpecialized<3301>(a);
    case 464:
        return RadicalInverseSpecialized<3307>(a);
    case 465:
        return RadicalInverseSpecialized<3313>(a);
    case 466:
        return RadicalInverseSpecialized<3319>(a);
    case 467:
        return RadicalInverseSpecialized<3323>(a);
    case 468:
        return RadicalInverseSpecialized<3329>(a);
    case 469:
        return RadicalInverseSpecialized<3331>(a);
    case 470:
        return RadicalInverseSpecialized<3343>(a);
    case 471:
        return RadicalInverseSpecialized<3347>(a);
    case 472:
        return RadicalInverseSpecialized<3359>(a);
    case 473:
        return RadicalInverseSpecialized<3361>(a);
    case 474:
        return RadicalInverseSpecialized<3371>(a);
    case 475:
        return RadicalInverseSpecialized<3373>(a);
    case 476:
        return RadicalInverseSpecialized<3389>(a);
    case 477:
        return RadicalInverseSpecialized<3391>(a);
    case 478:
        return RadicalInverseSpecialized<3407>(a);
    case 479:
        return RadicalInverseSpecialized<3413>(a);
    case 480:
        return RadicalInverseSpecialized<3433>(a);
    case 481:
        return RadicalInverseSpecialized<3449>(a);
    case 482:
        return RadicalInverseSpecialized<3457>(a);
    case 483:
        return RadicalInverseSpecialized<3461>(a);
    case 484:
        return RadicalInverseSpecialized<3463>(a);
    case 485:
        return RadicalInverseSpecialized<3467>(a);
    case 486:
        return RadicalInverseSpecialized<3469>(a);
    case 487:
        return RadicalInverseSpecialized<3491>(a);
    case 488:
        return RadicalInverseSpecialized<3499>(a);
    case 489:
        return RadicalInverseSpecialized<3511>(a);
    case 490:
        return RadicalInverseSpecialized<3517>(a);
    case 491:
        return RadicalInverseSpecialized<3527>(a);
    case 492:
        return RadicalInverseSpecialized<3529>(a);
    case 493:
        return RadicalInverseSpecialized<3533>(a);
    case 494:
        return RadicalInverseSpecialized<3539>(a);
    case 495:
        return RadicalInverseSpecialized<3541>(a);
    case 496:
        return RadicalInverseSpecialized<3547>(a);
    case 497:
        return RadicalInverseSpecialized<3557>(a);
    case 498:
        return RadicalInverseSpecialized<3559>(a);
    case 499:
        return RadicalInverseSpecialized<3571>(a);
    case 500:
        return RadicalInverseSpecialized<3581>(a);
    case 501:
        return RadicalInverseSpecialized<3583>(a);
    case 502:
        return RadicalInverseSpecialized<3593>(a);
    case 503:
        return RadicalInverseSpecialized<3607>(a);
    case 504:
        return RadicalInverseSpecialized<3613>(a);
    case 505:
        return RadicalInverseSpecialized<3617>(a);
    case 506:
        return RadicalInverseSpecialized<3623>(a);
    case 507:
        return RadicalInverseSpecialized<3631>(a);
    case 508:
        return RadicalInverseSpecialized<3637>(a);
    case 509:
        return RadicalInverseSpecialized<3643>(a);
    case 510:
        return RadicalInverseSpecialized<3659>(a);
    case 511:
        return RadicalInverseSpecialized<3671>(a);
    case 512:
        return RadicalInverseSpecialized<3673>(a);
    case 513:
        return RadicalInverseSpecialized<3677>(a);
    case 514:
        return RadicalInverseSpecialized<3691>(a);
    case 515:
        return RadicalInverseSpecialized<3697>(a);
    case 516:
        return RadicalInverseSpecialized<3701>(a);
    case 517:
        return RadicalInverseSpecialized<3709>(a);
    case 518:
        return RadicalInverseSpecialized<3719>(a);
    case 519:
        return RadicalInverseSpecialized<3727>(a);
    case 520:
        return RadicalInverseSpecialized<3733>(a);
    case 521:
        return RadicalInverseSpecialized<3739>(a);
    case 522:
        return RadicalInverseSpecialized<3761>(a);
    case 523:
        return RadicalInverseSpecialized<3767>(a);
    case 524:
        return RadicalInverseSpecialized<3769>(a);
    case 525:
        return RadicalInverseSpecialized<3779>(a);
    case 526:
        return RadicalInverseSpecialized<3793>(a);
    case 527:
        return RadicalInverseSpecialized<3797>(a);
    case 528:
        return RadicalInverseSpecialized<3803>(a);
    case 529:
        return RadicalInverseSpecialized<3821>(a);
    case 530:
        return RadicalInverseSpecialized<3823>(a);
    case 531:
        return RadicalInverseSpecialized<3833>(a);
    case 532:
        return RadicalInverseSpecialized<3847>(a);
    case 533:
        return RadicalInverseSpecialized<3851>(a);
    case 534:
        return RadicalInverseSpecialized<3853>(a);
    case 535:
        return RadicalInverseSpecialized<3863>(a);
    case 536:
        return RadicalInverseSpecialized<3877>(a);
    case 537:
        return RadicalInverseSpecialized<3881>(a);
    case 538:
        return RadicalInverseSpecialized<3889>(a);
    case 539:
        return RadicalInverseSpecialized<3907>(a);
    case 540:
        return RadicalInverseSpecialized<3911>(a);
    case 541:
        return RadicalInverseSpecialized<3917>(a);
    case 542:
        return RadicalInverseSpecialized<3919>(a);
    case 543:
        return RadicalInverseSpecialized<3923>(a);
    case 544:
        return RadicalInverseSpecialized<3929>(a);
    case 545:
        return RadicalInverseSpecialized<3931>(a);
    case 546:
        return RadicalInverseSpecialized<3943>(a);
    case 547:
        return RadicalInverseSpecialized<3947>(a);
    case 548:
        return RadicalInverseSpecialized<3967>(a);
    case 549:
        return RadicalInverseSpecialized<3989>(a);
    case 550:
        return RadicalInverseSpecialized<4001>(a);
    case 551:
        return RadicalInverseSpecialized<4003>(a);
    case 552:
        return RadicalInverseSpecialized<4007>(a);
    case 553:
        return RadicalInverseSpecialized<4013>(a);
    case 554:
        return RadicalInverseSpecialized<4019>(a);
    case 555:
        return RadicalInverseSpecialized<4021>(a);
    case 556:
        return RadicalInverseSpecialized<4027>(a);
    case 557:
        return RadicalInverseSpecialized<4049>(a);
    case 558:
        return RadicalInverseSpecialized<4051>(a);
    case 559:
        return RadicalInverseSpecialized<4057>(a);
    case 560:
        return RadicalInverseSpecialized<4073>(a);
    case 561:
        return RadicalInverseSpecialized<4079>(a);
    case 562:
        return RadicalInverseSpecialized<4091>(a);
    case 563:
        return RadicalInverseSpecialized<4093>(a);
    case 564:
        return RadicalInverseSpecialized<4099>(a);
    case 565:
        return RadicalInverseSpecialized<4111>(a);
    case 566:
        return RadicalInverseSpecialized<4127>(a);
    case 567:
        return RadicalInverseSpecialized<4129>(a);
    case 568:
        return RadicalInverseSpecialized<4133>(a);
    case 569:
        return RadicalInverseSpecialized<4139>(a);
    case 570:
        return RadicalInverseSpecialized<4153>(a);
    case 571:
        return RadicalInverseSpecialized<4157>(a);
    case 572:
        return RadicalInverseSpecialized<4159>(a);
    case 573:
        return RadicalInverseSpecialized<4177>(a);
    case 574:
        return RadicalInverseSpecialized<4201>(a);
    case 575:
        return RadicalInverseSpecialized<4211>(a);
    case 576:
        return RadicalInverseSpecialized<4217>(a);
    case 577:
        return RadicalInverseSpecialized<4219>(a);
    case 578:
        return RadicalInverseSpecialized<4229>(a);
    case 579:
        return RadicalInverseSpecialized<4231>(a);
    case 580:
        return RadicalInverseSpecialized<4241>(a);
    case 581:
        return RadicalInverseSpecialized<4243>(a);
    case 582:
        return RadicalInverseSpecialized<4253>(a);
    case 583:
        return RadicalInverseSpecialized<4259>(a);
    case 584:
        return RadicalInverseSpecialized<4261>(a);
    case 585:
        return RadicalInverseSpecialized<4271>(a);
    case 586:
        return RadicalInverseSpecialized<4273>(a);
    case 587:
        return RadicalInverseSpecialized<4283>(a);
    case 588:
        return RadicalInverseSpecialized<4289>(a);
    case 589:
        return RadicalInverseSpecialized<4297>(a);
    case 590:
        return RadicalInverseSpecialized<4327>(a);
    case 591:
        return RadicalInverseSpecialized<4337>(a);
    case 592:
        return RadicalInverseSpecialized<4339>(a);
    case 593:
        return RadicalInverseSpecialized<4349>(a);
    case 594:
        return RadicalInverseSpecialized<4357>(a);
    case 595:
        return RadicalInverseSpecialized<4363>(a);
    case 596:
        return RadicalInverseSpecialized<4373>(a);
    case 597:
        return RadicalInverseSpecialized<4391>(a);
    case 598:
        return RadicalInverseSpecialized<4397>(a);
    case 599:
        return RadicalInverseSpecialized<4409>(a);
    case 600:
        return RadicalInverseSpecialized<4421>(a);
    case 601:
        return RadicalInverseSpecialized<4423>(a);
    case 602:
        return RadicalInverseSpecialized<4441>(a);
    case 603:
        return RadicalInverseSpecialized<4447>(a);
    case 604:
        return RadicalInverseSpecialized<4451>(a);
    case 605:
        return RadicalInverseSpecialized<4457>(a);
    case 606:
        return RadicalInverseSpecialized<4463>(a);
    case 607:
        return RadicalInverseSpecialized<4481>(a);
    case 608:
        return RadicalInverseSpecialized<4483>(a);
    case 609:
        return RadicalInverseSpecialized<4493>(a);
    case 610:
        return RadicalInverseSpecialized<4507>(a);
    case 611:
        return RadicalInverseSpecialized<4513>(a);
    case 612:
        return RadicalInverseSpecialized<4517>(a);
    case 613:
        return RadicalInverseSpecialized<4519>(a);
    case 614:
        return RadicalInverseSpecialized<4523>(a);
    case 615:
        return RadicalInverseSpecialized<4547>(a);
    case 616:
        return RadicalInverseSpecialized<4549>(a);
    case 617:
        return RadicalInverseSpecialized<4561>(a);
    case 618:
        return RadicalInverseSpecialized<4567>(a);
    case 619:
        return RadicalInverseSpecialized<4583>(a);
    case 620:
        return RadicalInverseSpecialized<4591>(a);
    case 621:
        return RadicalInverseSpecialized<4597>(a);
    case 622:
        return RadicalInverseSpecialized<4603>(a);
    case 623:
        return RadicalInverseSpecialized<4621>(a);
    case 624:
        return RadicalInverseSpecialized<4637>(a);
    case 625:
        return RadicalInverseSpecialized<4639>(a);
    case 626:
        return RadicalInverseSpecialized<4643>(a);
    case 627:
        return RadicalInverseSpecialized<4649>(a);
    case 628:
        return RadicalInverseSpecialized<4651>(a);
    case 629:
        return RadicalInverseSpecialized<4657>(a);
    case 630:
        return RadicalInverseSpecialized<4663>(a);
    case 631:
        return RadicalInverseSpecialized<4673>(a);
    case 632:
        return RadicalInverseSpecialized<4679>(a);
    case 633:
        return RadicalInverseSpecialized<4691>(a);
    case 634:
        return RadicalInverseSpecialized<4703>(a);
    case 635:
        return RadicalInverseSpecialized<4721>(a);
    case 636:
        return RadicalInverseSpecialized<4723>(a);
    case 637:
        return RadicalInverseSpecialized<4729>(a);
    case 638:
        return RadicalInverseSpecialized<4733>(a);
    case 639:
        return RadicalInverseSpecialized<4751>(a);
    case 640:
        return RadicalInverseSpecialized<4759>(a);
    case 641:
        return RadicalInverseSpecialized<4783>(a);
    case 642:
        return RadicalInverseSpecialized<4787>(a);
    case 643:
        return RadicalInverseSpecialized<4789>(a);
    case 644:
        return RadicalInverseSpecialized<4793>(a);
    case 645:
        return RadicalInverseSpecialized<4799>(a);
    case 646:
        return RadicalInverseSpecialized<4801>(a);
    case 647:
        return RadicalInverseSpecialized<4813>(a);
    case 648:
        return RadicalInverseSpecialized<4817>(a);
    case 649:
        return RadicalInverseSpecialized<4831>(a);
    case 650:
        return RadicalInverseSpecialized<4861>(a);
    case 651:
        return RadicalInverseSpecialized<4871>(a);
    case 652:
        return RadicalInverseSpecialized<4877>(a);
    case 653:
        return RadicalInverseSpecialized<4889>(a);
    case 654:
        return RadicalInverseSpecialized<4903>(a);
    case 655:
        return RadicalInverseSpecialized<4909>(a);
    case 656:
        return RadicalInverseSpecialized<4919>(a);
    case 657:
        return RadicalInverseSpecialized<4931>(a);
    case 658:
        return RadicalInverseSpecialized<4933>(a);
    case 659:
        return RadicalInverseSpecialized<4937>(a);
    case 660:
        return RadicalInverseSpecialized<4943>(a);
    case 661:
        return RadicalInverseSpecialized<4951>(a);
    case 662:
        return RadicalInverseSpecialized<4957>(a);
    case 663:
        return RadicalInverseSpecialized<4967>(a);
    case 664:
        return RadicalInverseSpecialized<4969>(a);
    case 665:
        return RadicalInverseSpecialized<4973>(a);
    case 666:
        return RadicalInverseSpecialized<4987>(a);
    case 667:
        return RadicalInverseSpecialized<4993>(a);
    case 668:
        return RadicalInverseSpecialized<4999>(a);
    case 669:
        return RadicalInverseSpecialized<5003>(a);
    case 670:
        return RadicalInverseSpecialized<5009>(a);
    case 671:
        return RadicalInverseSpecialized<5011>(a);
    case 672:
        return RadicalInverseSpecialized<5021>(a);
    case 673:
        return RadicalInverseSpecialized<5023>(a);
    case 674:
        return RadicalInverseSpecialized<5039>(a);
    case 675:
        return RadicalInverseSpecialized<5051>(a);
    case 676:
        return RadicalInverseSpecialized<5059>(a);
    case 677:
        return RadicalInverseSpecialized<5077>(a);
    case 678:
        return RadicalInverseSpecialized<5081>(a);
    case 679:
        return RadicalInverseSpecialized<5087>(a);
    case 680:
        return RadicalInverseSpecialized<5099>(a);
    case 681:
        return RadicalInverseSpecialized<5101>(a);
    case 682:
        return RadicalInverseSpecialized<5107>(a);
    case 683:
        return RadicalInverseSpecialized<5113>(a);
    case 684:
        return RadicalInverseSpecialized<5119>(a);
    case 685:
        return RadicalInverseSpecialized<5147>(a);
    case 686:
        return RadicalInverseSpecialized<5153>(a);
    case 687:
        return RadicalInverseSpecialized<5167>(a);
    case 688:
        return RadicalInverseSpecialized<5171>(a);
    case 689:
        return RadicalInverseSpecialized<5179>(a);
    case 690:
        return RadicalInverseSpecialized<5189>(a);
    case 691:
        return RadicalInverseSpecialized<5197>(a);
    case 692:
        return RadicalInverseSpecialized<5209>(a);
    case 693:
        return RadicalInverseSpecialized<5227>(a);
    case 694:
        return RadicalInverseSpecialized<5231>(a);
    case 695:
        return RadicalInverseSpecialized<5233>(a);
    case 696:
        return RadicalInverseSpecialized<5237>(a);
    case 697:
        return RadicalInverseSpecialized<5261>(a);
    case 698:
        return RadicalInverseSpecialized<5273>(a);
    case 699:
        return RadicalInverseSpecialized<5279>(a);
    case 700:
        return RadicalInverseSpecialized<5281>(a);
    case 701:
        return RadicalInverseSpecialized<5297>(a);
    case 702:
        return RadicalInverseSpecialized<5303>(a);
    case 703:
        return RadicalInverseSpecialized<5309>(a);
    case 704:
        return RadicalInverseSpecialized<5323>(a);
    case 705:
        return RadicalInverseSpecialized<5333>(a);
    case 706:
        return RadicalInverseSpecialized<5347>(a);
    case 707:
        return RadicalInverseSpecialized<5351>(a);
    case 708:
        return RadicalInverseSpecialized<5381>(a);
    case 709:
        return RadicalInverseSpecialized<5387>(a);
    case 710:
        return RadicalInverseSpecialized<5393>(a);
    case 711:
        return RadicalInverseSpecialized<5399>(a);
    case 712:
        return RadicalInverseSpecialized<5407>(a);
    case 713:
        return RadicalInverseSpecialized<5413>(a);
    case 714:
        return RadicalInverseSpecialized<5417>(a);
    case 715:
        return RadicalInverseSpecialized<5419>(a);
    case 716:
        return RadicalInverseSpecialized<5431>(a);
    case 717:
        return RadicalInverseSpecialized<5437>(a);
    case 718:
        return RadicalInverseSpecialized<5441>(a);
    case 719:
        return RadicalInverseSpecialized<5443>(a);
    case 720:
        return RadicalInverseSpecialized<5449>(a);
    case 721:
        return RadicalInverseSpecialized<5471>(a);
    case 722:
        return RadicalInverseSpecialized<5477>(a);
    case 723:
        return RadicalInverseSpecialized<5479>(a);
    case 724:
        return RadicalInverseSpecialized<5483>(a);
    case 725:
        return RadicalInverseSpecialized<5501>(a);
    case 726:
        return RadicalInverseSpecialized<5503>(a);
    case 727:
        return RadicalInverseSpecialized<5507>(a);
    case 728:
        return RadicalInverseSpecialized<5519>(a);
    case 729:
        return RadicalInverseSpecialized<5521>(a);
    case 730:
        return RadicalInverseSpecialized<5527>(a);
    case 731:
        return RadicalInverseSpecialized<5531>(a);
    case 732:
        return RadicalInverseSpecialized<5557>(a);
    case 733:
        return RadicalInverseSpecialized<5563>(a);
    case 734:
        return RadicalInverseSpecialized<5569>(a);
    case 735:
        return RadicalInverseSpecialized<5573>(a);
    case 736:
        return RadicalInverseSpecialized<5581>(a);
    case 737:
        return RadicalInverseSpecialized<5591>(a);
    case 738:
        return RadicalInverseSpecialized<5623>(a);
    case 739:
        return RadicalInverseSpecialized<5639>(a);
    case 740:
        return RadicalInverseSpecialized<5641>(a);
    case 741:
        return RadicalInverseSpecialized<5647>(a);
    case 742:
        return RadicalInverseSpecialized<5651>(a);
    case 743:
        return RadicalInverseSpecialized<5653>(a);
    case 744:
        return RadicalInverseSpecialized<5657>(a);
    case 745:
        return RadicalInverseSpecialized<5659>(a);
    case 746:
        return RadicalInverseSpecialized<5669>(a);
    case 747:
        return RadicalInverseSpecialized<5683>(a);
    case 748:
        return RadicalInverseSpecialized<5689>(a);
    case 749:
        return RadicalInverseSpecialized<5693>(a);
    case 750:
        return RadicalInverseSpecialized<5701>(a);
    case 751:
        return RadicalInverseSpecialized<5711>(a);
    case 752:
        return RadicalInverseSpecialized<5717>(a);
    case 753:
        return RadicalInverseSpecialized<5737>(a);
    case 754:
        return RadicalInverseSpecialized<5741>(a);
    case 755:
        return RadicalInverseSpecialized<5743>(a);
    case 756:
        return RadicalInverseSpecialized<5749>(a);
    case 757:
        return RadicalInverseSpecialized<5779>(a);
    case 758:
        return RadicalInverseSpecialized<5783>(a);
    case 759:
        return RadicalInverseSpecialized<5791>(a);
    case 760:
        return RadicalInverseSpecialized<5801>(a);
    case 761:
        return RadicalInverseSpecialized<5807>(a);
    case 762:
        return RadicalInverseSpecialized<5813>(a);
    case 763:
        return RadicalInverseSpecialized<5821>(a);
    case 764:
        return RadicalInverseSpecialized<5827>(a);
    case 765:
        return RadicalInverseSpecialized<5839>(a);
    case 766:
        return RadicalInverseSpecialized<5843>(a);
    case 767:
        return RadicalInverseSpecialized<5849>(a);
    case 768:
        return RadicalInverseSpecialized<5851>(a);
    case 769:
        return RadicalInverseSpecialized<5857>(a);
    case 770:
        return RadicalInverseSpecialized<5861>(a);
    case 771:
        return RadicalInverseSpecialized<5867>(a);
    case 772:
        return RadicalInverseSpecialized<5869>(a);
    case 773:
        return RadicalInverseSpecialized<5879>(a);
    case 774:
        return RadicalInverseSpecialized<5881>(a);
    case 775:
        return RadicalInverseSpecialized<5897>(a);
    case 776:
        return RadicalInverseSpecialized<5903>(a);
    case 777:
        return RadicalInverseSpecialized<5923>(a);
    case 778:
        return RadicalInverseSpecialized<5927>(a);
    case 779:
        return RadicalInverseSpecialized<5939>(a);
    case 780:
        return RadicalInverseSpecialized<5953>(a);
    case 781:
        return RadicalInverseSpecialized<5981>(a);
    case 782:
        return RadicalInverseSpecialized<5987>(a);
    case 783:
        return RadicalInverseSpecialized<6007>(a);
    case 784:
        return RadicalInverseSpecialized<6011>(a);
    case 785:
        return RadicalInverseSpecialized<6029>(a);
    case 786:
        return RadicalInverseSpecialized<6037>(a);
    case 787:
        return RadicalInverseSpecialized<6043>(a);
    case 788:
        return RadicalInverseSpecialized<6047>(a);
    case 789:
        return RadicalInverseSpecialized<6053>(a);
    case 790:
        return RadicalInverseSpecialized<6067>(a);
    case 791:
        return RadicalInverseSpecialized<6073>(a);
    case 792:
        return RadicalInverseSpecialized<6079>(a);
    case 793:
        return RadicalInverseSpecialized<6089>(a);
    case 794:
        return RadicalInverseSpecialized<6091>(a);
    case 795:
        return RadicalInverseSpecialized<6101>(a);
    case 796:
        return RadicalInverseSpecialized<6113>(a);
    case 797:
        return RadicalInverseSpecialized<6121>(a);
    case 798:
        return RadicalInverseSpecialized<6131>(a);
    case 799:
        return RadicalInverseSpecialized<6133>(a);
    case 800:
        return RadicalInverseSpecialized<6143>(a);
    case 801:
        return RadicalInverseSpecialized<6151>(a);
    case 802:
        return RadicalInverseSpecialized<6163>(a);
    case 803:
        return RadicalInverseSpecialized<6173>(a);
    case 804:
        return RadicalInverseSpecialized<6197>(a);
    case 805:
        return RadicalInverseSpecialized<6199>(a);
    case 806:
        return RadicalInverseSpecialized<6203>(a);
    case 807:
        return RadicalInverseSpecialized<6211>(a);
    case 808:
        return RadicalInverseSpecialized<6217>(a);
    case 809:
        return RadicalInverseSpecialized<6221>(a);
    case 810:
        return RadicalInverseSpecialized<6229>(a);
    case 811:
        return RadicalInverseSpecialized<6247>(a);
    case 812:
        return RadicalInverseSpecialized<6257>(a);
    case 813:
        return RadicalInverseSpecialized<6263>(a);
    case 814:
        return RadicalInverseSpecialized<6269>(a);
    case 815:
        return RadicalInverseSpecialized<6271>(a);
    case 816:
        return RadicalInverseSpecialized<6277>(a);
    case 817:
        return RadicalInverseSpecialized<6287>(a);
    case 818:
        return RadicalInverseSpecialized<6299>(a);
    case 819:
        return RadicalInverseSpecialized<6301>(a);
    case 820:
        return RadicalInverseSpecialized<6311>(a);
    case 821:
        return RadicalInverseSpecialized<6317>(a);
    case 822:
        return RadicalInverseSpecialized<6323>(a);
    case 823:
        return RadicalInverseSpecialized<6329>(a);
    case 824:
        return RadicalInverseSpecialized<6337>(a);
    case 825:
        return RadicalInverseSpecialized<6343>(a);
    case 826:
        return RadicalInverseSpecialized<6353>(a);
    case 827:
        return RadicalInverseSpecialized<6359>(a);
    case 828:
        return RadicalInverseSpecialized<6361>(a);
    case 829:
        return RadicalInverseSpecialized<6367>(a);
    case 830:
        return RadicalInverseSpecialized<6373>(a);
    case 831:
        return RadicalInverseSpecialized<6379>(a);
    case 832:
        return RadicalInverseSpecialized<6389>(a);
    case 833:
        return RadicalInverseSpecialized<6397>(a);
    case 834:
        return RadicalInverseSpecialized<6421>(a);
    case 835:
        return RadicalInverseSpecialized<6427>(a);
    case 836:
        return RadicalInverseSpecialized<6449>(a);
    case 837:
        return RadicalInverseSpecialized<6451>(a);
    case 838:
        return RadicalInverseSpecialized<6469>(a);
    case 839:
        return RadicalInverseSpecialized<6473>(a);
    case 840:
        return RadicalInverseSpecialized<6481>(a);
    case 841:
        return RadicalInverseSpecialized<6491>(a);
    case 842:
        return RadicalInverseSpecialized<6521>(a);
    case 843:
        return RadicalInverseSpecialized<6529>(a);
    case 844:
        return RadicalInverseSpecialized<6547>(a);
    case 845:
        return RadicalInverseSpecialized<6551>(a);
    case 846:
        return RadicalInverseSpecialized<6553>(a);
    case 847:
        return RadicalInverseSpecialized<6563>(a);
    case 848:
        return RadicalInverseSpecialized<6569>(a);
    case 849:
        return RadicalInverseSpecialized<6571>(a);
    case 850:
        return RadicalInverseSpecialized<6577>(a);
    case 851:
        return RadicalInverseSpecialized<6581>(a);
    case 852:
        return RadicalInverseSpecialized<6599>(a);
    case 853:
        return RadicalInverseSpecialized<6607>(a);
    case 854:
        return RadicalInverseSpecialized<6619>(a);
    case 855:
        return RadicalInverseSpecialized<6637>(a);
    case 856:
        return RadicalInverseSpecialized<6653>(a);
    case 857:
        return RadicalInverseSpecialized<6659>(a);
    case 858:
        return RadicalInverseSpecialized<6661>(a);
    case 859:
        return RadicalInverseSpecialized<6673>(a);
    case 860:
        return RadicalInverseSpecialized<6679>(a);
    case 861:
        return RadicalInverseSpecialized<6689>(a);
    case 862:
        return RadicalInverseSpecialized<6691>(a);
    case 863:
        return RadicalInverseSpecialized<6701>(a);
    case 864:
        return RadicalInverseSpecialized<6703>(a);
    case 865:
        return RadicalInverseSpecialized<6709>(a);
    case 866:
        return RadicalInverseSpecialized<6719>(a);
    case 867:
        return RadicalInverseSpecialized<6733>(a);
    case 868:
        return RadicalInverseSpecialized<6737>(a);
    case 869:
        return RadicalInverseSpecialized<6761>(a);
    case 870:
        return RadicalInverseSpecialized<6763>(a);
    case 871:
        return RadicalInverseSpecialized<6779>(a);
    case 872:
        return RadicalInverseSpecialized<6781>(a);
    case 873:
        return RadicalInverseSpecialized<6791>(a);
    case 874:
        return RadicalInverseSpecialized<6793>(a);
    case 875:
        return RadicalInverseSpecialized<6803>(a);
    case 876:
        return RadicalInverseSpecialized<6823>(a);
    case 877:
        return RadicalInverseSpecialized<6827>(a);
    case 878:
        return RadicalInverseSpecialized<6829>(a);
    case 879:
        return RadicalInverseSpecialized<6833>(a);
    case 880:
        return RadicalInverseSpecialized<6841>(a);
    case 881:
        return RadicalInverseSpecialized<6857>(a);
    case 882:
        return RadicalInverseSpecialized<6863>(a);
    case 883:
        return RadicalInverseSpecialized<6869>(a);
    case 884:
        return RadicalInverseSpecialized<6871>(a);
    case 885:
        return RadicalInverseSpecialized<6883>(a);
    case 886:
        return RadicalInverseSpecialized<6899>(a);
    case 887:
        return RadicalInverseSpecialized<6907>(a);
    case 888:
        return RadicalInverseSpecialized<6911>(a);
    case 889:
        return RadicalInverseSpecialized<6917>(a);
    case 890:
        return RadicalInverseSpecialized<6947>(a);
    case 891:
        return RadicalInverseSpecialized<6949>(a);
    case 892:
        return RadicalInverseSpecialized<6959>(a);
    case 893:
        return RadicalInverseSpecialized<6961>(a);
    case 894:
        return RadicalInverseSpecialized<6967>(a);
    case 895:
        return RadicalInverseSpecialized<6971>(a);
    case 896:
        return RadicalInverseSpecialized<6977>(a);
    case 897:
        return RadicalInverseSpecialized<6983>(a);
    case 898:
        return RadicalInverseSpecialized<6991>(a);
    case 899:
        return RadicalInverseSpecialized<6997>(a);
    case 900:
        return RadicalInverseSpecialized<7001>(a);
    case 901:
        return RadicalInverseSpecialized<7013>(a);
    case 902:
        return RadicalInverseSpecialized<7019>(a);
    case 903:
        return RadicalInverseSpecialized<7027>(a);
    case 904:
        return RadicalInverseSpecialized<7039>(a);
    case 905:
        return RadicalInverseSpecialized<7043>(a);
    case 906:
        return RadicalInverseSpecialized<7057>(a);
    case 907:
        return RadicalInverseSpecialized<7069>(a);
    case 908:
        return RadicalInverseSpecialized<7079>(a);
    case 909:
        return RadicalInverseSpecialized<7103>(a);
    case 910:
        return RadicalInverseSpecialized<7109>(a);
    case 911:
        return RadicalInverseSpecialized<7121>(a);
    case 912:
        return RadicalInverseSpecialized<7127>(a);
    case 913:
        return RadicalInverseSpecialized<7129>(a);
    case 914:
        return RadicalInverseSpecialized<7151>(a);
    case 915:
        return RadicalInverseSpecialized<7159>(a);
    case 916:
        return RadicalInverseSpecialized<7177>(a);
    case 917:
        return RadicalInverseSpecialized<7187>(a);
    case 918:
        return RadicalInverseSpecialized<7193>(a);
    case 919:
        return RadicalInverseSpecialized<7207>(a);
    case 920:
        return RadicalInverseSpecialized<7211>(a);
    case 921:
        return RadicalInverseSpecialized<7213>(a);
    case 922:
        return RadicalInverseSpecialized<7219>(a);
    case 923:
        return RadicalInverseSpecialized<7229>(a);
    case 924:
        return RadicalInverseSpecialized<7237>(a);
    case 925:
        return RadicalInverseSpecialized<7243>(a);
    case 926:
        return RadicalInverseSpecialized<7247>(a);
    case 927:
        return RadicalInverseSpecialized<7253>(a);
    case 928:
        return RadicalInverseSpecialized<7283>(a);
    case 929:
        return RadicalInverseSpecialized<7297>(a);
    case 930:
        return RadicalInverseSpecialized<7307>(a);
    case 931:
        return RadicalInverseSpecialized<7309>(a);
    case 932:
        return RadicalInverseSpecialized<7321>(a);
    case 933:
        return RadicalInverseSpecialized<7331>(a);
    case 934:
        return RadicalInverseSpecialized<7333>(a);
    case 935:
        return RadicalInverseSpecialized<7349>(a);
    case 936:
        return RadicalInverseSpecialized<7351>(a);
    case 937:
        return RadicalInverseSpecialized<7369>(a);
    case 938:
        return RadicalInverseSpecialized<7393>(a);
    case 939:
        return RadicalInverseSpecialized<7411>(a);
    case 940:
        return RadicalInverseSpecialized<7417>(a);
    case 941:
        return RadicalInverseSpecialized<7433>(a);
    case 942:
        return RadicalInverseSpecialized<7451>(a);
    case 943:
        return RadicalInverseSpecialized<7457>(a);
    case 944:
        return RadicalInverseSpecialized<7459>(a);
    case 945:
        return RadicalInverseSpecialized<7477>(a);
    case 946:
        return RadicalInverseSpecialized<7481>(a);
    case 947:
        return RadicalInverseSpecialized<7487>(a);
    case 948:
        return RadicalInverseSpecialized<7489>(a);
    case 949:
        return RadicalInverseSpecialized<7499>(a);
    case 950:
        return RadicalInverseSpecialized<7507>(a);
    case 951:
        return RadicalInverseSpecialized<7517>(a);
    case 952:
        return RadicalInverseSpecialized<7523>(a);
    case 953:
        return RadicalInverseSpecialized<7529>(a);
    case 954:
        return RadicalInverseSpecialized<7537>(a);
    case 955:
        return RadicalInverseSpecialized<7541>(a);
    case 956:
        return RadicalInverseSpecialized<7547>(a);
    case 957:
        return RadicalInverseSpecialized<7549>(a);
    case 958:
        return RadicalInverseSpecialized<7559>(a);
    case 959:
        return RadicalInverseSpecialized<7561>(a);
    case 960:
        return RadicalInverseSpecialized<7573>(a);
    case 961:
        return RadicalInverseSpecialized<7577>(a);
    case 962:
        return RadicalInverseSpecialized<7583>(a);
    case 963:
        return RadicalInverseSpecialized<7589>(a);
    case 964:
        return RadicalInverseSpecialized<7591>(a);
    case 965:
        return RadicalInverseSpecialized<7603>(a);
    case 966:
        return RadicalInverseSpecialized<7607>(a);
    case 967:
        return RadicalInverseSpecialized<7621>(a);
    case 968:
        return RadicalInverseSpecialized<7639>(a);
    case 969:
        return RadicalInverseSpecialized<7643>(a);
    case 970:
        return RadicalInverseSpecialized<7649>(a);
    case 971:
        return RadicalInverseSpecialized<7669>(a);
    case 972:
        return RadicalInverseSpecialized<7673>(a);
    case 973:
        return RadicalInverseSpecialized<7681>(a);
    case 974:
        return RadicalInverseSpecialized<7687>(a);
    case 975:
        return RadicalInverseSpecialized<7691>(a);
    case 976:
        return RadicalInverseSpecialized<7699>(a);
    case 977:
        return RadicalInverseSpecialized<7703>(a);
    case 978:
        return RadicalInverseSpecialized<7717>(a);
    case 979:
        return RadicalInverseSpecialized<7723>(a);
    case 980:
        return RadicalInverseSpecialized<7727>(a);
    case 981:
        return RadicalInverseSpecialized<7741>(a);
    case 982:
        return RadicalInverseSpecialized<7753>(a);
    case 983:
        return RadicalInverseSpecialized<7757>(a);
    case 984:
        return RadicalInverseSpecialized<7759>(a);
    case 985:
        return RadicalInverseSpecialized<7789>(a);
    case 986:
        return RadicalInverseSpecialized<7793>(a);
    case 987:
        return RadicalInverseSpecialized<7817>(a);
    case 988:
        return RadicalInverseSpecialized<7823>(a);
    case 989:
        return RadicalInverseSpecialized<7829>(a);
    case 990:
        return RadicalInverseSpecialized<7841>(a);
    case 991:
        return RadicalInverseSpecialized<7853>(a);
    case 992:
        return RadicalInverseSpecialized<7867>(a);
    case 993:
        return RadicalInverseSpecialized<7873>(a);
    case 994:
        return RadicalInverseSpecialized<7877>(a);
    case 995:
        return RadicalInverseSpecialized<7879>(a);
    case 996:
        return RadicalInverseSpecialized<7883>(a);
    case 997:
        return RadicalInverseSpecialized<7901>(a);
    case 998:
        return RadicalInverseSpecialized<7907>(a);
    case 999:
        return RadicalInverseSpecialized<7919>(a);
    case 1000:
        return RadicalInverseSpecialized<7927>(a);
    case 1001:
        return RadicalInverseSpecialized<7933>(a);
    case 1002:
        return RadicalInverseSpecialized<7937>(a);
    case 1003:
        return RadicalInverseSpecialized<7949>(a);
    case 1004:
        return RadicalInverseSpecialized<7951>(a);
    case 1005:
        return RadicalInverseSpecialized<7963>(a);
    case 1006:
        return RadicalInverseSpecialized<7993>(a);
    case 1007:
        return RadicalInverseSpecialized<8009>(a);
    case 1008:
        return RadicalInverseSpecialized<8011>(a);
    case 1009:
        return RadicalInverseSpecialized<8017>(a);
    case 1010:
        return RadicalInverseSpecialized<8039>(a);
    case 1011:
        return RadicalInverseSpecialized<8053>(a);
    case 1012:
        return RadicalInverseSpecialized<8059>(a);
    case 1013:
        return RadicalInverseSpecialized<8069>(a);
    case 1014:
        return RadicalInverseSpecialized<8081>(a);
    case 1015:
        return RadicalInverseSpecialized<8087>(a);
    case 1016:
        return RadicalInverseSpecialized<8089>(a);
    case 1017:
        return RadicalInverseSpecialized<8093>(a);
    case 1018:
        return RadicalInverseSpecialized<8101>(a);
    case 1019:
        return RadicalInverseSpecialized<8111>(a);
    case 1020:
        return RadicalInverseSpecialized<8117>(a);
    case 1021:
        return RadicalInverseSpecialized<8123>(a);
    case 1022:
        return RadicalInverseSpecialized<8147>(a);
    case 1023:
        return RadicalInverseSpecialized<8161>(a);
    default:
        // LOG(FATAL) << StringPrintf("Base %d is >= 1024, the limit of RadicalInverse",
        //                            baseIndex);
        return 0;
    }
}

__both__
Float ScrambledRadicalInverse(int baseIndex, uint64_t a, const uint16_t *perm) {
    switch (baseIndex) {
    case 0:
        return ScrambledRadicalInverseSpecialized<2>(perm, a);
    case 1:
        return ScrambledRadicalInverseSpecialized<3>(perm, a);
    case 2:
        return ScrambledRadicalInverseSpecialized<5>(perm, a);
    case 3:
        return ScrambledRadicalInverseSpecialized<7>(perm, a);
    // Remainder of cases for _ScrambledRadicalInverse()_
    case 4:
        return ScrambledRadicalInverseSpecialized<11>(perm, a);
    case 5:
        return ScrambledRadicalInverseSpecialized<13>(perm, a);
    case 6:
        return ScrambledRadicalInverseSpecialized<17>(perm, a);
    case 7:
        return ScrambledRadicalInverseSpecialized<19>(perm, a);
    case 8:
        return ScrambledRadicalInverseSpecialized<23>(perm, a);
    case 9:
        return ScrambledRadicalInverseSpecialized<29>(perm, a);
    case 10:
        return ScrambledRadicalInverseSpecialized<31>(perm, a);
    case 11:
        return ScrambledRadicalInverseSpecialized<37>(perm, a);
    case 12:
        return ScrambledRadicalInverseSpecialized<41>(perm, a);
    case 13:
        return ScrambledRadicalInverseSpecialized<43>(perm, a);
    case 14:
        return ScrambledRadicalInverseSpecialized<47>(perm, a);
    case 15:
        return ScrambledRadicalInverseSpecialized<53>(perm, a);
    case 16:
        return ScrambledRadicalInverseSpecialized<59>(perm, a);
    case 17:
        return ScrambledRadicalInverseSpecialized<61>(perm, a);
    case 18:
        return ScrambledRadicalInverseSpecialized<67>(perm, a);
    case 19:
        return ScrambledRadicalInverseSpecialized<71>(perm, a);
    case 20:
        return ScrambledRadicalInverseSpecialized<73>(perm, a);
    case 21:
        return ScrambledRadicalInverseSpecialized<79>(perm, a);
    case 22:
        return ScrambledRadicalInverseSpecialized<83>(perm, a);
    case 23:
        return ScrambledRadicalInverseSpecialized<89>(perm, a);
    case 24:
        return ScrambledRadicalInverseSpecialized<97>(perm, a);
    case 25:
        return ScrambledRadicalInverseSpecialized<101>(perm, a);
    case 26:
        return ScrambledRadicalInverseSpecialized<103>(perm, a);
    case 27:
        return ScrambledRadicalInverseSpecialized<107>(perm, a);
    case 28:
        return ScrambledRadicalInverseSpecialized<109>(perm, a);
    case 29:
        return ScrambledRadicalInverseSpecialized<113>(perm, a);
    case 30:
        return ScrambledRadicalInverseSpecialized<127>(perm, a);
    case 31:
        return ScrambledRadicalInverseSpecialized<131>(perm, a);
    case 32:
        return ScrambledRadicalInverseSpecialized<137>(perm, a);
    case 33:
        return ScrambledRadicalInverseSpecialized<139>(perm, a);
    case 34:
        return ScrambledRadicalInverseSpecialized<149>(perm, a);
    case 35:
        return ScrambledRadicalInverseSpecialized<151>(perm, a);
    case 36:
        return ScrambledRadicalInverseSpecialized<157>(perm, a);
    case 37:
        return ScrambledRadicalInverseSpecialized<163>(perm, a);
    case 38:
        return ScrambledRadicalInverseSpecialized<167>(perm, a);
    case 39:
        return ScrambledRadicalInverseSpecialized<173>(perm, a);
    case 40:
        return ScrambledRadicalInverseSpecialized<179>(perm, a);
    case 41:
        return ScrambledRadicalInverseSpecialized<181>(perm, a);
    case 42:
        return ScrambledRadicalInverseSpecialized<191>(perm, a);
    case 43:
        return ScrambledRadicalInverseSpecialized<193>(perm, a);
    case 44:
        return ScrambledRadicalInverseSpecialized<197>(perm, a);
    case 45:
        return ScrambledRadicalInverseSpecialized<199>(perm, a);
    case 46:
        return ScrambledRadicalInverseSpecialized<211>(perm, a);
    case 47:
        return ScrambledRadicalInverseSpecialized<223>(perm, a);
    case 48:
        return ScrambledRadicalInverseSpecialized<227>(perm, a);
    case 49:
        return ScrambledRadicalInverseSpecialized<229>(perm, a);
    case 50:
        return ScrambledRadicalInverseSpecialized<233>(perm, a);
    case 51:
        return ScrambledRadicalInverseSpecialized<239>(perm, a);
    case 52:
        return ScrambledRadicalInverseSpecialized<241>(perm, a);
    case 53:
        return ScrambledRadicalInverseSpecialized<251>(perm, a);
    case 54:
        return ScrambledRadicalInverseSpecialized<257>(perm, a);
    case 55:
        return ScrambledRadicalInverseSpecialized<263>(perm, a);
    case 56:
        return ScrambledRadicalInverseSpecialized<269>(perm, a);
    case 57:
        return ScrambledRadicalInverseSpecialized<271>(perm, a);
    case 58:
        return ScrambledRadicalInverseSpecialized<277>(perm, a);
    case 59:
        return ScrambledRadicalInverseSpecialized<281>(perm, a);
    case 60:
        return ScrambledRadicalInverseSpecialized<283>(perm, a);
    case 61:
        return ScrambledRadicalInverseSpecialized<293>(perm, a);
    case 62:
        return ScrambledRadicalInverseSpecialized<307>(perm, a);
    case 63:
        return ScrambledRadicalInverseSpecialized<311>(perm, a);
    case 64:
        return ScrambledRadicalInverseSpecialized<313>(perm, a);
    case 65:
        return ScrambledRadicalInverseSpecialized<317>(perm, a);
    case 66:
        return ScrambledRadicalInverseSpecialized<331>(perm, a);
    case 67:
        return ScrambledRadicalInverseSpecialized<337>(perm, a);
    case 68:
        return ScrambledRadicalInverseSpecialized<347>(perm, a);
    case 69:
        return ScrambledRadicalInverseSpecialized<349>(perm, a);
    case 70:
        return ScrambledRadicalInverseSpecialized<353>(perm, a);
    case 71:
        return ScrambledRadicalInverseSpecialized<359>(perm, a);
    case 72:
        return ScrambledRadicalInverseSpecialized<367>(perm, a);
    case 73:
        return ScrambledRadicalInverseSpecialized<373>(perm, a);
    case 74:
        return ScrambledRadicalInverseSpecialized<379>(perm, a);
    case 75:
        return ScrambledRadicalInverseSpecialized<383>(perm, a);
    case 76:
        return ScrambledRadicalInverseSpecialized<389>(perm, a);
    case 77:
        return ScrambledRadicalInverseSpecialized<397>(perm, a);
    case 78:
        return ScrambledRadicalInverseSpecialized<401>(perm, a);
    case 79:
        return ScrambledRadicalInverseSpecialized<409>(perm, a);
    case 80:
        return ScrambledRadicalInverseSpecialized<419>(perm, a);
    case 81:
        return ScrambledRadicalInverseSpecialized<421>(perm, a);
    case 82:
        return ScrambledRadicalInverseSpecialized<431>(perm, a);
    case 83:
        return ScrambledRadicalInverseSpecialized<433>(perm, a);
    case 84:
        return ScrambledRadicalInverseSpecialized<439>(perm, a);
    case 85:
        return ScrambledRadicalInverseSpecialized<443>(perm, a);
    case 86:
        return ScrambledRadicalInverseSpecialized<449>(perm, a);
    case 87:
        return ScrambledRadicalInverseSpecialized<457>(perm, a);
    case 88:
        return ScrambledRadicalInverseSpecialized<461>(perm, a);
    case 89:
        return ScrambledRadicalInverseSpecialized<463>(perm, a);
    case 90:
        return ScrambledRadicalInverseSpecialized<467>(perm, a);
    case 91:
        return ScrambledRadicalInverseSpecialized<479>(perm, a);
    case 92:
        return ScrambledRadicalInverseSpecialized<487>(perm, a);
    case 93:
        return ScrambledRadicalInverseSpecialized<491>(perm, a);
    case 94:
        return ScrambledRadicalInverseSpecialized<499>(perm, a);
    case 95:
        return ScrambledRadicalInverseSpecialized<503>(perm, a);
    case 96:
        return ScrambledRadicalInverseSpecialized<509>(perm, a);
    case 97:
        return ScrambledRadicalInverseSpecialized<521>(perm, a);
    case 98:
        return ScrambledRadicalInverseSpecialized<523>(perm, a);
    case 99:
        return ScrambledRadicalInverseSpecialized<541>(perm, a);
    case 100:
        return ScrambledRadicalInverseSpecialized<547>(perm, a);
    case 101:
        return ScrambledRadicalInverseSpecialized<557>(perm, a);
    case 102:
        return ScrambledRadicalInverseSpecialized<563>(perm, a);
    case 103:
        return ScrambledRadicalInverseSpecialized<569>(perm, a);
    case 104:
        return ScrambledRadicalInverseSpecialized<571>(perm, a);
    case 105:
        return ScrambledRadicalInverseSpecialized<577>(perm, a);
    case 106:
        return ScrambledRadicalInverseSpecialized<587>(perm, a);
    case 107:
        return ScrambledRadicalInverseSpecialized<593>(perm, a);
    case 108:
        return ScrambledRadicalInverseSpecialized<599>(perm, a);
    case 109:
        return ScrambledRadicalInverseSpecialized<601>(perm, a);
    case 110:
        return ScrambledRadicalInverseSpecialized<607>(perm, a);
    case 111:
        return ScrambledRadicalInverseSpecialized<613>(perm, a);
    case 112:
        return ScrambledRadicalInverseSpecialized<617>(perm, a);
    case 113:
        return ScrambledRadicalInverseSpecialized<619>(perm, a);
    case 114:
        return ScrambledRadicalInverseSpecialized<631>(perm, a);
    case 115:
        return ScrambledRadicalInverseSpecialized<641>(perm, a);
    case 116:
        return ScrambledRadicalInverseSpecialized<643>(perm, a);
    case 117:
        return ScrambledRadicalInverseSpecialized<647>(perm, a);
    case 118:
        return ScrambledRadicalInverseSpecialized<653>(perm, a);
    case 119:
        return ScrambledRadicalInverseSpecialized<659>(perm, a);
    case 120:
        return ScrambledRadicalInverseSpecialized<661>(perm, a);
    case 121:
        return ScrambledRadicalInverseSpecialized<673>(perm, a);
    case 122:
        return ScrambledRadicalInverseSpecialized<677>(perm, a);
    case 123:
        return ScrambledRadicalInverseSpecialized<683>(perm, a);
    case 124:
        return ScrambledRadicalInverseSpecialized<691>(perm, a);
    case 125:
        return ScrambledRadicalInverseSpecialized<701>(perm, a);
    case 126:
        return ScrambledRadicalInverseSpecialized<709>(perm, a);
    case 127:
        return ScrambledRadicalInverseSpecialized<719>(perm, a);
    case 128:
        return ScrambledRadicalInverseSpecialized<727>(perm, a);
    case 129:
        return ScrambledRadicalInverseSpecialized<733>(perm, a);
    case 130:
        return ScrambledRadicalInverseSpecialized<739>(perm, a);
    case 131:
        return ScrambledRadicalInverseSpecialized<743>(perm, a);
    case 132:
        return ScrambledRadicalInverseSpecialized<751>(perm, a);
    case 133:
        return ScrambledRadicalInverseSpecialized<757>(perm, a);
    case 134:
        return ScrambledRadicalInverseSpecialized<761>(perm, a);
    case 135:
        return ScrambledRadicalInverseSpecialized<769>(perm, a);
    case 136:
        return ScrambledRadicalInverseSpecialized<773>(perm, a);
    case 137:
        return ScrambledRadicalInverseSpecialized<787>(perm, a);
    case 138:
        return ScrambledRadicalInverseSpecialized<797>(perm, a);
    case 139:
        return ScrambledRadicalInverseSpecialized<809>(perm, a);
    case 140:
        return ScrambledRadicalInverseSpecialized<811>(perm, a);
    case 141:
        return ScrambledRadicalInverseSpecialized<821>(perm, a);
    case 142:
        return ScrambledRadicalInverseSpecialized<823>(perm, a);
    case 143:
        return ScrambledRadicalInverseSpecialized<827>(perm, a);
    case 144:
        return ScrambledRadicalInverseSpecialized<829>(perm, a);
    case 145:
        return ScrambledRadicalInverseSpecialized<839>(perm, a);
    case 146:
        return ScrambledRadicalInverseSpecialized<853>(perm, a);
    case 147:
        return ScrambledRadicalInverseSpecialized<857>(perm, a);
    case 148:
        return ScrambledRadicalInverseSpecialized<859>(perm, a);
    case 149:
        return ScrambledRadicalInverseSpecialized<863>(perm, a);
    case 150:
        return ScrambledRadicalInverseSpecialized<877>(perm, a);
    case 151:
        return ScrambledRadicalInverseSpecialized<881>(perm, a);
    case 152:
        return ScrambledRadicalInverseSpecialized<883>(perm, a);
    case 153:
        return ScrambledRadicalInverseSpecialized<887>(perm, a);
    case 154:
        return ScrambledRadicalInverseSpecialized<907>(perm, a);
    case 155:
        return ScrambledRadicalInverseSpecialized<911>(perm, a);
    case 156:
        return ScrambledRadicalInverseSpecialized<919>(perm, a);
    case 157:
        return ScrambledRadicalInverseSpecialized<929>(perm, a);
    case 158:
        return ScrambledRadicalInverseSpecialized<937>(perm, a);
    case 159:
        return ScrambledRadicalInverseSpecialized<941>(perm, a);
    case 160:
        return ScrambledRadicalInverseSpecialized<947>(perm, a);
    case 161:
        return ScrambledRadicalInverseSpecialized<953>(perm, a);
    case 162:
        return ScrambledRadicalInverseSpecialized<967>(perm, a);
    case 163:
        return ScrambledRadicalInverseSpecialized<971>(perm, a);
    case 164:
        return ScrambledRadicalInverseSpecialized<977>(perm, a);
    case 165:
        return ScrambledRadicalInverseSpecialized<983>(perm, a);
    case 166:
        return ScrambledRadicalInverseSpecialized<991>(perm, a);
    case 167:
        return ScrambledRadicalInverseSpecialized<997>(perm, a);
    case 168:
        return ScrambledRadicalInverseSpecialized<1009>(perm, a);
    case 169:
        return ScrambledRadicalInverseSpecialized<1013>(perm, a);
    case 170:
        return ScrambledRadicalInverseSpecialized<1019>(perm, a);
    case 171:
        return ScrambledRadicalInverseSpecialized<1021>(perm, a);
    case 172:
        return ScrambledRadicalInverseSpecialized<1031>(perm, a);
    case 173:
        return ScrambledRadicalInverseSpecialized<1033>(perm, a);
    case 174:
        return ScrambledRadicalInverseSpecialized<1039>(perm, a);
    case 175:
        return ScrambledRadicalInverseSpecialized<1049>(perm, a);
    case 176:
        return ScrambledRadicalInverseSpecialized<1051>(perm, a);
    case 177:
        return ScrambledRadicalInverseSpecialized<1061>(perm, a);
    case 178:
        return ScrambledRadicalInverseSpecialized<1063>(perm, a);
    case 179:
        return ScrambledRadicalInverseSpecialized<1069>(perm, a);
    case 180:
        return ScrambledRadicalInverseSpecialized<1087>(perm, a);
    case 181:
        return ScrambledRadicalInverseSpecialized<1091>(perm, a);
    case 182:
        return ScrambledRadicalInverseSpecialized<1093>(perm, a);
    case 183:
        return ScrambledRadicalInverseSpecialized<1097>(perm, a);
    case 184:
        return ScrambledRadicalInverseSpecialized<1103>(perm, a);
    case 185:
        return ScrambledRadicalInverseSpecialized<1109>(perm, a);
    case 186:
        return ScrambledRadicalInverseSpecialized<1117>(perm, a);
    case 187:
        return ScrambledRadicalInverseSpecialized<1123>(perm, a);
    case 188:
        return ScrambledRadicalInverseSpecialized<1129>(perm, a);
    case 189:
        return ScrambledRadicalInverseSpecialized<1151>(perm, a);
    case 190:
        return ScrambledRadicalInverseSpecialized<1153>(perm, a);
    case 191:
        return ScrambledRadicalInverseSpecialized<1163>(perm, a);
    case 192:
        return ScrambledRadicalInverseSpecialized<1171>(perm, a);
    case 193:
        return ScrambledRadicalInverseSpecialized<1181>(perm, a);
    case 194:
        return ScrambledRadicalInverseSpecialized<1187>(perm, a);
    case 195:
        return ScrambledRadicalInverseSpecialized<1193>(perm, a);
    case 196:
        return ScrambledRadicalInverseSpecialized<1201>(perm, a);
    case 197:
        return ScrambledRadicalInverseSpecialized<1213>(perm, a);
    case 198:
        return ScrambledRadicalInverseSpecialized<1217>(perm, a);
    case 199:
        return ScrambledRadicalInverseSpecialized<1223>(perm, a);
    case 200:
        return ScrambledRadicalInverseSpecialized<1229>(perm, a);
    case 201:
        return ScrambledRadicalInverseSpecialized<1231>(perm, a);
    case 202:
        return ScrambledRadicalInverseSpecialized<1237>(perm, a);
    case 203:
        return ScrambledRadicalInverseSpecialized<1249>(perm, a);
    case 204:
        return ScrambledRadicalInverseSpecialized<1259>(perm, a);
    case 205:
        return ScrambledRadicalInverseSpecialized<1277>(perm, a);
    case 206:
        return ScrambledRadicalInverseSpecialized<1279>(perm, a);
    case 207:
        return ScrambledRadicalInverseSpecialized<1283>(perm, a);
    case 208:
        return ScrambledRadicalInverseSpecialized<1289>(perm, a);
    case 209:
        return ScrambledRadicalInverseSpecialized<1291>(perm, a);
    case 210:
        return ScrambledRadicalInverseSpecialized<1297>(perm, a);
    case 211:
        return ScrambledRadicalInverseSpecialized<1301>(perm, a);
    case 212:
        return ScrambledRadicalInverseSpecialized<1303>(perm, a);
    case 213:
        return ScrambledRadicalInverseSpecialized<1307>(perm, a);
    case 214:
        return ScrambledRadicalInverseSpecialized<1319>(perm, a);
    case 215:
        return ScrambledRadicalInverseSpecialized<1321>(perm, a);
    case 216:
        return ScrambledRadicalInverseSpecialized<1327>(perm, a);
    case 217:
        return ScrambledRadicalInverseSpecialized<1361>(perm, a);
    case 218:
        return ScrambledRadicalInverseSpecialized<1367>(perm, a);
    case 219:
        return ScrambledRadicalInverseSpecialized<1373>(perm, a);
    case 220:
        return ScrambledRadicalInverseSpecialized<1381>(perm, a);
    case 221:
        return ScrambledRadicalInverseSpecialized<1399>(perm, a);
    case 222:
        return ScrambledRadicalInverseSpecialized<1409>(perm, a);
    case 223:
        return ScrambledRadicalInverseSpecialized<1423>(perm, a);
    case 224:
        return ScrambledRadicalInverseSpecialized<1427>(perm, a);
    case 225:
        return ScrambledRadicalInverseSpecialized<1429>(perm, a);
    case 226:
        return ScrambledRadicalInverseSpecialized<1433>(perm, a);
    case 227:
        return ScrambledRadicalInverseSpecialized<1439>(perm, a);
    case 228:
        return ScrambledRadicalInverseSpecialized<1447>(perm, a);
    case 229:
        return ScrambledRadicalInverseSpecialized<1451>(perm, a);
    case 230:
        return ScrambledRadicalInverseSpecialized<1453>(perm, a);
    case 231:
        return ScrambledRadicalInverseSpecialized<1459>(perm, a);
    case 232:
        return ScrambledRadicalInverseSpecialized<1471>(perm, a);
    case 233:
        return ScrambledRadicalInverseSpecialized<1481>(perm, a);
    case 234:
        return ScrambledRadicalInverseSpecialized<1483>(perm, a);
    case 235:
        return ScrambledRadicalInverseSpecialized<1487>(perm, a);
    case 236:
        return ScrambledRadicalInverseSpecialized<1489>(perm, a);
    case 237:
        return ScrambledRadicalInverseSpecialized<1493>(perm, a);
    case 238:
        return ScrambledRadicalInverseSpecialized<1499>(perm, a);
    case 239:
        return ScrambledRadicalInverseSpecialized<1511>(perm, a);
    case 240:
        return ScrambledRadicalInverseSpecialized<1523>(perm, a);
    case 241:
        return ScrambledRadicalInverseSpecialized<1531>(perm, a);
    case 242:
        return ScrambledRadicalInverseSpecialized<1543>(perm, a);
    case 243:
        return ScrambledRadicalInverseSpecialized<1549>(perm, a);
    case 244:
        return ScrambledRadicalInverseSpecialized<1553>(perm, a);
    case 245:
        return ScrambledRadicalInverseSpecialized<1559>(perm, a);
    case 246:
        return ScrambledRadicalInverseSpecialized<1567>(perm, a);
    case 247:
        return ScrambledRadicalInverseSpecialized<1571>(perm, a);
    case 248:
        return ScrambledRadicalInverseSpecialized<1579>(perm, a);
    case 249:
        return ScrambledRadicalInverseSpecialized<1583>(perm, a);
    case 250:
        return ScrambledRadicalInverseSpecialized<1597>(perm, a);
    case 251:
        return ScrambledRadicalInverseSpecialized<1601>(perm, a);
    case 252:
        return ScrambledRadicalInverseSpecialized<1607>(perm, a);
    case 253:
        return ScrambledRadicalInverseSpecialized<1609>(perm, a);
    case 254:
        return ScrambledRadicalInverseSpecialized<1613>(perm, a);
    case 255:
        return ScrambledRadicalInverseSpecialized<1619>(perm, a);
    case 256:
        return ScrambledRadicalInverseSpecialized<1621>(perm, a);
    case 257:
        return ScrambledRadicalInverseSpecialized<1627>(perm, a);
    case 258:
        return ScrambledRadicalInverseSpecialized<1637>(perm, a);
    case 259:
        return ScrambledRadicalInverseSpecialized<1657>(perm, a);
    case 260:
        return ScrambledRadicalInverseSpecialized<1663>(perm, a);
    case 261:
        return ScrambledRadicalInverseSpecialized<1667>(perm, a);
    case 262:
        return ScrambledRadicalInverseSpecialized<1669>(perm, a);
    case 263:
        return ScrambledRadicalInverseSpecialized<1693>(perm, a);
    case 264:
        return ScrambledRadicalInverseSpecialized<1697>(perm, a);
    case 265:
        return ScrambledRadicalInverseSpecialized<1699>(perm, a);
    case 266:
        return ScrambledRadicalInverseSpecialized<1709>(perm, a);
    case 267:
        return ScrambledRadicalInverseSpecialized<1721>(perm, a);
    case 268:
        return ScrambledRadicalInverseSpecialized<1723>(perm, a);
    case 269:
        return ScrambledRadicalInverseSpecialized<1733>(perm, a);
    case 270:
        return ScrambledRadicalInverseSpecialized<1741>(perm, a);
    case 271:
        return ScrambledRadicalInverseSpecialized<1747>(perm, a);
    case 272:
        return ScrambledRadicalInverseSpecialized<1753>(perm, a);
    case 273:
        return ScrambledRadicalInverseSpecialized<1759>(perm, a);
    case 274:
        return ScrambledRadicalInverseSpecialized<1777>(perm, a);
    case 275:
        return ScrambledRadicalInverseSpecialized<1783>(perm, a);
    case 276:
        return ScrambledRadicalInverseSpecialized<1787>(perm, a);
    case 277:
        return ScrambledRadicalInverseSpecialized<1789>(perm, a);
    case 278:
        return ScrambledRadicalInverseSpecialized<1801>(perm, a);
    case 279:
        return ScrambledRadicalInverseSpecialized<1811>(perm, a);
    case 280:
        return ScrambledRadicalInverseSpecialized<1823>(perm, a);
    case 281:
        return ScrambledRadicalInverseSpecialized<1831>(perm, a);
    case 282:
        return ScrambledRadicalInverseSpecialized<1847>(perm, a);
    case 283:
        return ScrambledRadicalInverseSpecialized<1861>(perm, a);
    case 284:
        return ScrambledRadicalInverseSpecialized<1867>(perm, a);
    case 285:
        return ScrambledRadicalInverseSpecialized<1871>(perm, a);
    case 286:
        return ScrambledRadicalInverseSpecialized<1873>(perm, a);
    case 287:
        return ScrambledRadicalInverseSpecialized<1877>(perm, a);
    case 288:
        return ScrambledRadicalInverseSpecialized<1879>(perm, a);
    case 289:
        return ScrambledRadicalInverseSpecialized<1889>(perm, a);
    case 290:
        return ScrambledRadicalInverseSpecialized<1901>(perm, a);
    case 291:
        return ScrambledRadicalInverseSpecialized<1907>(perm, a);
    case 292:
        return ScrambledRadicalInverseSpecialized<1913>(perm, a);
    case 293:
        return ScrambledRadicalInverseSpecialized<1931>(perm, a);
    case 294:
        return ScrambledRadicalInverseSpecialized<1933>(perm, a);
    case 295:
        return ScrambledRadicalInverseSpecialized<1949>(perm, a);
    case 296:
        return ScrambledRadicalInverseSpecialized<1951>(perm, a);
    case 297:
        return ScrambledRadicalInverseSpecialized<1973>(perm, a);
    case 298:
        return ScrambledRadicalInverseSpecialized<1979>(perm, a);
    case 299:
        return ScrambledRadicalInverseSpecialized<1987>(perm, a);
    case 300:
        return ScrambledRadicalInverseSpecialized<1993>(perm, a);
    case 301:
        return ScrambledRadicalInverseSpecialized<1997>(perm, a);
    case 302:
        return ScrambledRadicalInverseSpecialized<1999>(perm, a);
    case 303:
        return ScrambledRadicalInverseSpecialized<2003>(perm, a);
    case 304:
        return ScrambledRadicalInverseSpecialized<2011>(perm, a);
    case 305:
        return ScrambledRadicalInverseSpecialized<2017>(perm, a);
    case 306:
        return ScrambledRadicalInverseSpecialized<2027>(perm, a);
    case 307:
        return ScrambledRadicalInverseSpecialized<2029>(perm, a);
    case 308:
        return ScrambledRadicalInverseSpecialized<2039>(perm, a);
    case 309:
        return ScrambledRadicalInverseSpecialized<2053>(perm, a);
    case 310:
        return ScrambledRadicalInverseSpecialized<2063>(perm, a);
    case 311:
        return ScrambledRadicalInverseSpecialized<2069>(perm, a);
    case 312:
        return ScrambledRadicalInverseSpecialized<2081>(perm, a);
    case 313:
        return ScrambledRadicalInverseSpecialized<2083>(perm, a);
    case 314:
        return ScrambledRadicalInverseSpecialized<2087>(perm, a);
    case 315:
        return ScrambledRadicalInverseSpecialized<2089>(perm, a);
    case 316:
        return ScrambledRadicalInverseSpecialized<2099>(perm, a);
    case 317:
        return ScrambledRadicalInverseSpecialized<2111>(perm, a);
    case 318:
        return ScrambledRadicalInverseSpecialized<2113>(perm, a);
    case 319:
        return ScrambledRadicalInverseSpecialized<2129>(perm, a);
    case 320:
        return ScrambledRadicalInverseSpecialized<2131>(perm, a);
    case 321:
        return ScrambledRadicalInverseSpecialized<2137>(perm, a);
    case 322:
        return ScrambledRadicalInverseSpecialized<2141>(perm, a);
    case 323:
        return ScrambledRadicalInverseSpecialized<2143>(perm, a);
    case 324:
        return ScrambledRadicalInverseSpecialized<2153>(perm, a);
    case 325:
        return ScrambledRadicalInverseSpecialized<2161>(perm, a);
    case 326:
        return ScrambledRadicalInverseSpecialized<2179>(perm, a);
    case 327:
        return ScrambledRadicalInverseSpecialized<2203>(perm, a);
    case 328:
        return ScrambledRadicalInverseSpecialized<2207>(perm, a);
    case 329:
        return ScrambledRadicalInverseSpecialized<2213>(perm, a);
    case 330:
        return ScrambledRadicalInverseSpecialized<2221>(perm, a);
    case 331:
        return ScrambledRadicalInverseSpecialized<2237>(perm, a);
    case 332:
        return ScrambledRadicalInverseSpecialized<2239>(perm, a);
    case 333:
        return ScrambledRadicalInverseSpecialized<2243>(perm, a);
    case 334:
        return ScrambledRadicalInverseSpecialized<2251>(perm, a);
    case 335:
        return ScrambledRadicalInverseSpecialized<2267>(perm, a);
    case 336:
        return ScrambledRadicalInverseSpecialized<2269>(perm, a);
    case 337:
        return ScrambledRadicalInverseSpecialized<2273>(perm, a);
    case 338:
        return ScrambledRadicalInverseSpecialized<2281>(perm, a);
    case 339:
        return ScrambledRadicalInverseSpecialized<2287>(perm, a);
    case 340:
        return ScrambledRadicalInverseSpecialized<2293>(perm, a);
    case 341:
        return ScrambledRadicalInverseSpecialized<2297>(perm, a);
    case 342:
        return ScrambledRadicalInverseSpecialized<2309>(perm, a);
    case 343:
        return ScrambledRadicalInverseSpecialized<2311>(perm, a);
    case 344:
        return ScrambledRadicalInverseSpecialized<2333>(perm, a);
    case 345:
        return ScrambledRadicalInverseSpecialized<2339>(perm, a);
    case 346:
        return ScrambledRadicalInverseSpecialized<2341>(perm, a);
    case 347:
        return ScrambledRadicalInverseSpecialized<2347>(perm, a);
    case 348:
        return ScrambledRadicalInverseSpecialized<2351>(perm, a);
    case 349:
        return ScrambledRadicalInverseSpecialized<2357>(perm, a);
    case 350:
        return ScrambledRadicalInverseSpecialized<2371>(perm, a);
    case 351:
        return ScrambledRadicalInverseSpecialized<2377>(perm, a);
    case 352:
        return ScrambledRadicalInverseSpecialized<2381>(perm, a);
    case 353:
        return ScrambledRadicalInverseSpecialized<2383>(perm, a);
    case 354:
        return ScrambledRadicalInverseSpecialized<2389>(perm, a);
    case 355:
        return ScrambledRadicalInverseSpecialized<2393>(perm, a);
    case 356:
        return ScrambledRadicalInverseSpecialized<2399>(perm, a);
    case 357:
        return ScrambledRadicalInverseSpecialized<2411>(perm, a);
    case 358:
        return ScrambledRadicalInverseSpecialized<2417>(perm, a);
    case 359:
        return ScrambledRadicalInverseSpecialized<2423>(perm, a);
    case 360:
        return ScrambledRadicalInverseSpecialized<2437>(perm, a);
    case 361:
        return ScrambledRadicalInverseSpecialized<2441>(perm, a);
    case 362:
        return ScrambledRadicalInverseSpecialized<2447>(perm, a);
    case 363:
        return ScrambledRadicalInverseSpecialized<2459>(perm, a);
    case 364:
        return ScrambledRadicalInverseSpecialized<2467>(perm, a);
    case 365:
        return ScrambledRadicalInverseSpecialized<2473>(perm, a);
    case 366:
        return ScrambledRadicalInverseSpecialized<2477>(perm, a);
    case 367:
        return ScrambledRadicalInverseSpecialized<2503>(perm, a);
    case 368:
        return ScrambledRadicalInverseSpecialized<2521>(perm, a);
    case 369:
        return ScrambledRadicalInverseSpecialized<2531>(perm, a);
    case 370:
        return ScrambledRadicalInverseSpecialized<2539>(perm, a);
    case 371:
        return ScrambledRadicalInverseSpecialized<2543>(perm, a);
    case 372:
        return ScrambledRadicalInverseSpecialized<2549>(perm, a);
    case 373:
        return ScrambledRadicalInverseSpecialized<2551>(perm, a);
    case 374:
        return ScrambledRadicalInverseSpecialized<2557>(perm, a);
    case 375:
        return ScrambledRadicalInverseSpecialized<2579>(perm, a);
    case 376:
        return ScrambledRadicalInverseSpecialized<2591>(perm, a);
    case 377:
        return ScrambledRadicalInverseSpecialized<2593>(perm, a);
    case 378:
        return ScrambledRadicalInverseSpecialized<2609>(perm, a);
    case 379:
        return ScrambledRadicalInverseSpecialized<2617>(perm, a);
    case 380:
        return ScrambledRadicalInverseSpecialized<2621>(perm, a);
    case 381:
        return ScrambledRadicalInverseSpecialized<2633>(perm, a);
    case 382:
        return ScrambledRadicalInverseSpecialized<2647>(perm, a);
    case 383:
        return ScrambledRadicalInverseSpecialized<2657>(perm, a);
    case 384:
        return ScrambledRadicalInverseSpecialized<2659>(perm, a);
    case 385:
        return ScrambledRadicalInverseSpecialized<2663>(perm, a);
    case 386:
        return ScrambledRadicalInverseSpecialized<2671>(perm, a);
    case 387:
        return ScrambledRadicalInverseSpecialized<2677>(perm, a);
    case 388:
        return ScrambledRadicalInverseSpecialized<2683>(perm, a);
    case 389:
        return ScrambledRadicalInverseSpecialized<2687>(perm, a);
    case 390:
        return ScrambledRadicalInverseSpecialized<2689>(perm, a);
    case 391:
        return ScrambledRadicalInverseSpecialized<2693>(perm, a);
    case 392:
        return ScrambledRadicalInverseSpecialized<2699>(perm, a);
    case 393:
        return ScrambledRadicalInverseSpecialized<2707>(perm, a);
    case 394:
        return ScrambledRadicalInverseSpecialized<2711>(perm, a);
    case 395:
        return ScrambledRadicalInverseSpecialized<2713>(perm, a);
    case 396:
        return ScrambledRadicalInverseSpecialized<2719>(perm, a);
    case 397:
        return ScrambledRadicalInverseSpecialized<2729>(perm, a);
    case 398:
        return ScrambledRadicalInverseSpecialized<2731>(perm, a);
    case 399:
        return ScrambledRadicalInverseSpecialized<2741>(perm, a);
    case 400:
        return ScrambledRadicalInverseSpecialized<2749>(perm, a);
    case 401:
        return ScrambledRadicalInverseSpecialized<2753>(perm, a);
    case 402:
        return ScrambledRadicalInverseSpecialized<2767>(perm, a);
    case 403:
        return ScrambledRadicalInverseSpecialized<2777>(perm, a);
    case 404:
        return ScrambledRadicalInverseSpecialized<2789>(perm, a);
    case 405:
        return ScrambledRadicalInverseSpecialized<2791>(perm, a);
    case 406:
        return ScrambledRadicalInverseSpecialized<2797>(perm, a);
    case 407:
        return ScrambledRadicalInverseSpecialized<2801>(perm, a);
    case 408:
        return ScrambledRadicalInverseSpecialized<2803>(perm, a);
    case 409:
        return ScrambledRadicalInverseSpecialized<2819>(perm, a);
    case 410:
        return ScrambledRadicalInverseSpecialized<2833>(perm, a);
    case 411:
        return ScrambledRadicalInverseSpecialized<2837>(perm, a);
    case 412:
        return ScrambledRadicalInverseSpecialized<2843>(perm, a);
    case 413:
        return ScrambledRadicalInverseSpecialized<2851>(perm, a);
    case 414:
        return ScrambledRadicalInverseSpecialized<2857>(perm, a);
    case 415:
        return ScrambledRadicalInverseSpecialized<2861>(perm, a);
    case 416:
        return ScrambledRadicalInverseSpecialized<2879>(perm, a);
    case 417:
        return ScrambledRadicalInverseSpecialized<2887>(perm, a);
    case 418:
        return ScrambledRadicalInverseSpecialized<2897>(perm, a);
    case 419:
        return ScrambledRadicalInverseSpecialized<2903>(perm, a);
    case 420:
        return ScrambledRadicalInverseSpecialized<2909>(perm, a);
    case 421:
        return ScrambledRadicalInverseSpecialized<2917>(perm, a);
    case 422:
        return ScrambledRadicalInverseSpecialized<2927>(perm, a);
    case 423:
        return ScrambledRadicalInverseSpecialized<2939>(perm, a);
    case 424:
        return ScrambledRadicalInverseSpecialized<2953>(perm, a);
    case 425:
        return ScrambledRadicalInverseSpecialized<2957>(perm, a);
    case 426:
        return ScrambledRadicalInverseSpecialized<2963>(perm, a);
    case 427:
        return ScrambledRadicalInverseSpecialized<2969>(perm, a);
    case 428:
        return ScrambledRadicalInverseSpecialized<2971>(perm, a);
    case 429:
        return ScrambledRadicalInverseSpecialized<2999>(perm, a);
    case 430:
        return ScrambledRadicalInverseSpecialized<3001>(perm, a);
    case 431:
        return ScrambledRadicalInverseSpecialized<3011>(perm, a);
    case 432:
        return ScrambledRadicalInverseSpecialized<3019>(perm, a);
    case 433:
        return ScrambledRadicalInverseSpecialized<3023>(perm, a);
    case 434:
        return ScrambledRadicalInverseSpecialized<3037>(perm, a);
    case 435:
        return ScrambledRadicalInverseSpecialized<3041>(perm, a);
    case 436:
        return ScrambledRadicalInverseSpecialized<3049>(perm, a);
    case 437:
        return ScrambledRadicalInverseSpecialized<3061>(perm, a);
    case 438:
        return ScrambledRadicalInverseSpecialized<3067>(perm, a);
    case 439:
        return ScrambledRadicalInverseSpecialized<3079>(perm, a);
    case 440:
        return ScrambledRadicalInverseSpecialized<3083>(perm, a);
    case 441:
        return ScrambledRadicalInverseSpecialized<3089>(perm, a);
    case 442:
        return ScrambledRadicalInverseSpecialized<3109>(perm, a);
    case 443:
        return ScrambledRadicalInverseSpecialized<3119>(perm, a);
    case 444:
        return ScrambledRadicalInverseSpecialized<3121>(perm, a);
    case 445:
        return ScrambledRadicalInverseSpecialized<3137>(perm, a);
    case 446:
        return ScrambledRadicalInverseSpecialized<3163>(perm, a);
    case 447:
        return ScrambledRadicalInverseSpecialized<3167>(perm, a);
    case 448:
        return ScrambledRadicalInverseSpecialized<3169>(perm, a);
    case 449:
        return ScrambledRadicalInverseSpecialized<3181>(perm, a);
    case 450:
        return ScrambledRadicalInverseSpecialized<3187>(perm, a);
    case 451:
        return ScrambledRadicalInverseSpecialized<3191>(perm, a);
    case 452:
        return ScrambledRadicalInverseSpecialized<3203>(perm, a);
    case 453:
        return ScrambledRadicalInverseSpecialized<3209>(perm, a);
    case 454:
        return ScrambledRadicalInverseSpecialized<3217>(perm, a);
    case 455:
        return ScrambledRadicalInverseSpecialized<3221>(perm, a);
    case 456:
        return ScrambledRadicalInverseSpecialized<3229>(perm, a);
    case 457:
        return ScrambledRadicalInverseSpecialized<3251>(perm, a);
    case 458:
        return ScrambledRadicalInverseSpecialized<3253>(perm, a);
    case 459:
        return ScrambledRadicalInverseSpecialized<3257>(perm, a);
    case 460:
        return ScrambledRadicalInverseSpecialized<3259>(perm, a);
    case 461:
        return ScrambledRadicalInverseSpecialized<3271>(perm, a);
    case 462:
        return ScrambledRadicalInverseSpecialized<3299>(perm, a);
    case 463:
        return ScrambledRadicalInverseSpecialized<3301>(perm, a);
    case 464:
        return ScrambledRadicalInverseSpecialized<3307>(perm, a);
    case 465:
        return ScrambledRadicalInverseSpecialized<3313>(perm, a);
    case 466:
        return ScrambledRadicalInverseSpecialized<3319>(perm, a);
    case 467:
        return ScrambledRadicalInverseSpecialized<3323>(perm, a);
    case 468:
        return ScrambledRadicalInverseSpecialized<3329>(perm, a);
    case 469:
        return ScrambledRadicalInverseSpecialized<3331>(perm, a);
    case 470:
        return ScrambledRadicalInverseSpecialized<3343>(perm, a);
    case 471:
        return ScrambledRadicalInverseSpecialized<3347>(perm, a);
    case 472:
        return ScrambledRadicalInverseSpecialized<3359>(perm, a);
    case 473:
        return ScrambledRadicalInverseSpecialized<3361>(perm, a);
    case 474:
        return ScrambledRadicalInverseSpecialized<3371>(perm, a);
    case 475:
        return ScrambledRadicalInverseSpecialized<3373>(perm, a);
    case 476:
        return ScrambledRadicalInverseSpecialized<3389>(perm, a);
    case 477:
        return ScrambledRadicalInverseSpecialized<3391>(perm, a);
    case 478:
        return ScrambledRadicalInverseSpecialized<3407>(perm, a);
    case 479:
        return ScrambledRadicalInverseSpecialized<3413>(perm, a);
    case 480:
        return ScrambledRadicalInverseSpecialized<3433>(perm, a);
    case 481:
        return ScrambledRadicalInverseSpecialized<3449>(perm, a);
    case 482:
        return ScrambledRadicalInverseSpecialized<3457>(perm, a);
    case 483:
        return ScrambledRadicalInverseSpecialized<3461>(perm, a);
    case 484:
        return ScrambledRadicalInverseSpecialized<3463>(perm, a);
    case 485:
        return ScrambledRadicalInverseSpecialized<3467>(perm, a);
    case 486:
        return ScrambledRadicalInverseSpecialized<3469>(perm, a);
    case 487:
        return ScrambledRadicalInverseSpecialized<3491>(perm, a);
    case 488:
        return ScrambledRadicalInverseSpecialized<3499>(perm, a);
    case 489:
        return ScrambledRadicalInverseSpecialized<3511>(perm, a);
    case 490:
        return ScrambledRadicalInverseSpecialized<3517>(perm, a);
    case 491:
        return ScrambledRadicalInverseSpecialized<3527>(perm, a);
    case 492:
        return ScrambledRadicalInverseSpecialized<3529>(perm, a);
    case 493:
        return ScrambledRadicalInverseSpecialized<3533>(perm, a);
    case 494:
        return ScrambledRadicalInverseSpecialized<3539>(perm, a);
    case 495:
        return ScrambledRadicalInverseSpecialized<3541>(perm, a);
    case 496:
        return ScrambledRadicalInverseSpecialized<3547>(perm, a);
    case 497:
        return ScrambledRadicalInverseSpecialized<3557>(perm, a);
    case 498:
        return ScrambledRadicalInverseSpecialized<3559>(perm, a);
    case 499:
        return ScrambledRadicalInverseSpecialized<3571>(perm, a);
    case 500:
        return ScrambledRadicalInverseSpecialized<3581>(perm, a);
    case 501:
        return ScrambledRadicalInverseSpecialized<3583>(perm, a);
    case 502:
        return ScrambledRadicalInverseSpecialized<3593>(perm, a);
    case 503:
        return ScrambledRadicalInverseSpecialized<3607>(perm, a);
    case 504:
        return ScrambledRadicalInverseSpecialized<3613>(perm, a);
    case 505:
        return ScrambledRadicalInverseSpecialized<3617>(perm, a);
    case 506:
        return ScrambledRadicalInverseSpecialized<3623>(perm, a);
    case 507:
        return ScrambledRadicalInverseSpecialized<3631>(perm, a);
    case 508:
        return ScrambledRadicalInverseSpecialized<3637>(perm, a);
    case 509:
        return ScrambledRadicalInverseSpecialized<3643>(perm, a);
    case 510:
        return ScrambledRadicalInverseSpecialized<3659>(perm, a);
    case 511:
        return ScrambledRadicalInverseSpecialized<3671>(perm, a);
    case 512:
        return ScrambledRadicalInverseSpecialized<3673>(perm, a);
    case 513:
        return ScrambledRadicalInverseSpecialized<3677>(perm, a);
    case 514:
        return ScrambledRadicalInverseSpecialized<3691>(perm, a);
    case 515:
        return ScrambledRadicalInverseSpecialized<3697>(perm, a);
    case 516:
        return ScrambledRadicalInverseSpecialized<3701>(perm, a);
    case 517:
        return ScrambledRadicalInverseSpecialized<3709>(perm, a);
    case 518:
        return ScrambledRadicalInverseSpecialized<3719>(perm, a);
    case 519:
        return ScrambledRadicalInverseSpecialized<3727>(perm, a);
    case 520:
        return ScrambledRadicalInverseSpecialized<3733>(perm, a);
    case 521:
        return ScrambledRadicalInverseSpecialized<3739>(perm, a);
    case 522:
        return ScrambledRadicalInverseSpecialized<3761>(perm, a);
    case 523:
        return ScrambledRadicalInverseSpecialized<3767>(perm, a);
    case 524:
        return ScrambledRadicalInverseSpecialized<3769>(perm, a);
    case 525:
        return ScrambledRadicalInverseSpecialized<3779>(perm, a);
    case 526:
        return ScrambledRadicalInverseSpecialized<3793>(perm, a);
    case 527:
        return ScrambledRadicalInverseSpecialized<3797>(perm, a);
    case 528:
        return ScrambledRadicalInverseSpecialized<3803>(perm, a);
    case 529:
        return ScrambledRadicalInverseSpecialized<3821>(perm, a);
    case 530:
        return ScrambledRadicalInverseSpecialized<3823>(perm, a);
    case 531:
        return ScrambledRadicalInverseSpecialized<3833>(perm, a);
    case 532:
        return ScrambledRadicalInverseSpecialized<3847>(perm, a);
    case 533:
        return ScrambledRadicalInverseSpecialized<3851>(perm, a);
    case 534:
        return ScrambledRadicalInverseSpecialized<3853>(perm, a);
    case 535:
        return ScrambledRadicalInverseSpecialized<3863>(perm, a);
    case 536:
        return ScrambledRadicalInverseSpecialized<3877>(perm, a);
    case 537:
        return ScrambledRadicalInverseSpecialized<3881>(perm, a);
    case 538:
        return ScrambledRadicalInverseSpecialized<3889>(perm, a);
    case 539:
        return ScrambledRadicalInverseSpecialized<3907>(perm, a);
    case 540:
        return ScrambledRadicalInverseSpecialized<3911>(perm, a);
    case 541:
        return ScrambledRadicalInverseSpecialized<3917>(perm, a);
    case 542:
        return ScrambledRadicalInverseSpecialized<3919>(perm, a);
    case 543:
        return ScrambledRadicalInverseSpecialized<3923>(perm, a);
    case 544:
        return ScrambledRadicalInverseSpecialized<3929>(perm, a);
    case 545:
        return ScrambledRadicalInverseSpecialized<3931>(perm, a);
    case 546:
        return ScrambledRadicalInverseSpecialized<3943>(perm, a);
    case 547:
        return ScrambledRadicalInverseSpecialized<3947>(perm, a);
    case 548:
        return ScrambledRadicalInverseSpecialized<3967>(perm, a);
    case 549:
        return ScrambledRadicalInverseSpecialized<3989>(perm, a);
    case 550:
        return ScrambledRadicalInverseSpecialized<4001>(perm, a);
    case 551:
        return ScrambledRadicalInverseSpecialized<4003>(perm, a);
    case 552:
        return ScrambledRadicalInverseSpecialized<4007>(perm, a);
    case 553:
        return ScrambledRadicalInverseSpecialized<4013>(perm, a);
    case 554:
        return ScrambledRadicalInverseSpecialized<4019>(perm, a);
    case 555:
        return ScrambledRadicalInverseSpecialized<4021>(perm, a);
    case 556:
        return ScrambledRadicalInverseSpecialized<4027>(perm, a);
    case 557:
        return ScrambledRadicalInverseSpecialized<4049>(perm, a);
    case 558:
        return ScrambledRadicalInverseSpecialized<4051>(perm, a);
    case 559:
        return ScrambledRadicalInverseSpecialized<4057>(perm, a);
    case 560:
        return ScrambledRadicalInverseSpecialized<4073>(perm, a);
    case 561:
        return ScrambledRadicalInverseSpecialized<4079>(perm, a);
    case 562:
        return ScrambledRadicalInverseSpecialized<4091>(perm, a);
    case 563:
        return ScrambledRadicalInverseSpecialized<4093>(perm, a);
    case 564:
        return ScrambledRadicalInverseSpecialized<4099>(perm, a);
    case 565:
        return ScrambledRadicalInverseSpecialized<4111>(perm, a);
    case 566:
        return ScrambledRadicalInverseSpecialized<4127>(perm, a);
    case 567:
        return ScrambledRadicalInverseSpecialized<4129>(perm, a);
    case 568:
        return ScrambledRadicalInverseSpecialized<4133>(perm, a);
    case 569:
        return ScrambledRadicalInverseSpecialized<4139>(perm, a);
    case 570:
        return ScrambledRadicalInverseSpecialized<4153>(perm, a);
    case 571:
        return ScrambledRadicalInverseSpecialized<4157>(perm, a);
    case 572:
        return ScrambledRadicalInverseSpecialized<4159>(perm, a);
    case 573:
        return ScrambledRadicalInverseSpecialized<4177>(perm, a);
    case 574:
        return ScrambledRadicalInverseSpecialized<4201>(perm, a);
    case 575:
        return ScrambledRadicalInverseSpecialized<4211>(perm, a);
    case 576:
        return ScrambledRadicalInverseSpecialized<4217>(perm, a);
    case 577:
        return ScrambledRadicalInverseSpecialized<4219>(perm, a);
    case 578:
        return ScrambledRadicalInverseSpecialized<4229>(perm, a);
    case 579:
        return ScrambledRadicalInverseSpecialized<4231>(perm, a);
    case 580:
        return ScrambledRadicalInverseSpecialized<4241>(perm, a);
    case 581:
        return ScrambledRadicalInverseSpecialized<4243>(perm, a);
    case 582:
        return ScrambledRadicalInverseSpecialized<4253>(perm, a);
    case 583:
        return ScrambledRadicalInverseSpecialized<4259>(perm, a);
    case 584:
        return ScrambledRadicalInverseSpecialized<4261>(perm, a);
    case 585:
        return ScrambledRadicalInverseSpecialized<4271>(perm, a);
    case 586:
        return ScrambledRadicalInverseSpecialized<4273>(perm, a);
    case 587:
        return ScrambledRadicalInverseSpecialized<4283>(perm, a);
    case 588:
        return ScrambledRadicalInverseSpecialized<4289>(perm, a);
    case 589:
        return ScrambledRadicalInverseSpecialized<4297>(perm, a);
    case 590:
        return ScrambledRadicalInverseSpecialized<4327>(perm, a);
    case 591:
        return ScrambledRadicalInverseSpecialized<4337>(perm, a);
    case 592:
        return ScrambledRadicalInverseSpecialized<4339>(perm, a);
    case 593:
        return ScrambledRadicalInverseSpecialized<4349>(perm, a);
    case 594:
        return ScrambledRadicalInverseSpecialized<4357>(perm, a);
    case 595:
        return ScrambledRadicalInverseSpecialized<4363>(perm, a);
    case 596:
        return ScrambledRadicalInverseSpecialized<4373>(perm, a);
    case 597:
        return ScrambledRadicalInverseSpecialized<4391>(perm, a);
    case 598:
        return ScrambledRadicalInverseSpecialized<4397>(perm, a);
    case 599:
        return ScrambledRadicalInverseSpecialized<4409>(perm, a);
    case 600:
        return ScrambledRadicalInverseSpecialized<4421>(perm, a);
    case 601:
        return ScrambledRadicalInverseSpecialized<4423>(perm, a);
    case 602:
        return ScrambledRadicalInverseSpecialized<4441>(perm, a);
    case 603:
        return ScrambledRadicalInverseSpecialized<4447>(perm, a);
    case 604:
        return ScrambledRadicalInverseSpecialized<4451>(perm, a);
    case 605:
        return ScrambledRadicalInverseSpecialized<4457>(perm, a);
    case 606:
        return ScrambledRadicalInverseSpecialized<4463>(perm, a);
    case 607:
        return ScrambledRadicalInverseSpecialized<4481>(perm, a);
    case 608:
        return ScrambledRadicalInverseSpecialized<4483>(perm, a);
    case 609:
        return ScrambledRadicalInverseSpecialized<4493>(perm, a);
    case 610:
        return ScrambledRadicalInverseSpecialized<4507>(perm, a);
    case 611:
        return ScrambledRadicalInverseSpecialized<4513>(perm, a);
    case 612:
        return ScrambledRadicalInverseSpecialized<4517>(perm, a);
    case 613:
        return ScrambledRadicalInverseSpecialized<4519>(perm, a);
    case 614:
        return ScrambledRadicalInverseSpecialized<4523>(perm, a);
    case 615:
        return ScrambledRadicalInverseSpecialized<4547>(perm, a);
    case 616:
        return ScrambledRadicalInverseSpecialized<4549>(perm, a);
    case 617:
        return ScrambledRadicalInverseSpecialized<4561>(perm, a);
    case 618:
        return ScrambledRadicalInverseSpecialized<4567>(perm, a);
    case 619:
        return ScrambledRadicalInverseSpecialized<4583>(perm, a);
    case 620:
        return ScrambledRadicalInverseSpecialized<4591>(perm, a);
    case 621:
        return ScrambledRadicalInverseSpecialized<4597>(perm, a);
    case 622:
        return ScrambledRadicalInverseSpecialized<4603>(perm, a);
    case 623:
        return ScrambledRadicalInverseSpecialized<4621>(perm, a);
    case 624:
        return ScrambledRadicalInverseSpecialized<4637>(perm, a);
    case 625:
        return ScrambledRadicalInverseSpecialized<4639>(perm, a);
    case 626:
        return ScrambledRadicalInverseSpecialized<4643>(perm, a);
    case 627:
        return ScrambledRadicalInverseSpecialized<4649>(perm, a);
    case 628:
        return ScrambledRadicalInverseSpecialized<4651>(perm, a);
    case 629:
        return ScrambledRadicalInverseSpecialized<4657>(perm, a);
    case 630:
        return ScrambledRadicalInverseSpecialized<4663>(perm, a);
    case 631:
        return ScrambledRadicalInverseSpecialized<4673>(perm, a);
    case 632:
        return ScrambledRadicalInverseSpecialized<4679>(perm, a);
    case 633:
        return ScrambledRadicalInverseSpecialized<4691>(perm, a);
    case 634:
        return ScrambledRadicalInverseSpecialized<4703>(perm, a);
    case 635:
        return ScrambledRadicalInverseSpecialized<4721>(perm, a);
    case 636:
        return ScrambledRadicalInverseSpecialized<4723>(perm, a);
    case 637:
        return ScrambledRadicalInverseSpecialized<4729>(perm, a);
    case 638:
        return ScrambledRadicalInverseSpecialized<4733>(perm, a);
    case 639:
        return ScrambledRadicalInverseSpecialized<4751>(perm, a);
    case 640:
        return ScrambledRadicalInverseSpecialized<4759>(perm, a);
    case 641:
        return ScrambledRadicalInverseSpecialized<4783>(perm, a);
    case 642:
        return ScrambledRadicalInverseSpecialized<4787>(perm, a);
    case 643:
        return ScrambledRadicalInverseSpecialized<4789>(perm, a);
    case 644:
        return ScrambledRadicalInverseSpecialized<4793>(perm, a);
    case 645:
        return ScrambledRadicalInverseSpecialized<4799>(perm, a);
    case 646:
        return ScrambledRadicalInverseSpecialized<4801>(perm, a);
    case 647:
        return ScrambledRadicalInverseSpecialized<4813>(perm, a);
    case 648:
        return ScrambledRadicalInverseSpecialized<4817>(perm, a);
    case 649:
        return ScrambledRadicalInverseSpecialized<4831>(perm, a);
    case 650:
        return ScrambledRadicalInverseSpecialized<4861>(perm, a);
    case 651:
        return ScrambledRadicalInverseSpecialized<4871>(perm, a);
    case 652:
        return ScrambledRadicalInverseSpecialized<4877>(perm, a);
    case 653:
        return ScrambledRadicalInverseSpecialized<4889>(perm, a);
    case 654:
        return ScrambledRadicalInverseSpecialized<4903>(perm, a);
    case 655:
        return ScrambledRadicalInverseSpecialized<4909>(perm, a);
    case 656:
        return ScrambledRadicalInverseSpecialized<4919>(perm, a);
    case 657:
        return ScrambledRadicalInverseSpecialized<4931>(perm, a);
    case 658:
        return ScrambledRadicalInverseSpecialized<4933>(perm, a);
    case 659:
        return ScrambledRadicalInverseSpecialized<4937>(perm, a);
    case 660:
        return ScrambledRadicalInverseSpecialized<4943>(perm, a);
    case 661:
        return ScrambledRadicalInverseSpecialized<4951>(perm, a);
    case 662:
        return ScrambledRadicalInverseSpecialized<4957>(perm, a);
    case 663:
        return ScrambledRadicalInverseSpecialized<4967>(perm, a);
    case 664:
        return ScrambledRadicalInverseSpecialized<4969>(perm, a);
    case 665:
        return ScrambledRadicalInverseSpecialized<4973>(perm, a);
    case 666:
        return ScrambledRadicalInverseSpecialized<4987>(perm, a);
    case 667:
        return ScrambledRadicalInverseSpecialized<4993>(perm, a);
    case 668:
        return ScrambledRadicalInverseSpecialized<4999>(perm, a);
    case 669:
        return ScrambledRadicalInverseSpecialized<5003>(perm, a);
    case 670:
        return ScrambledRadicalInverseSpecialized<5009>(perm, a);
    case 671:
        return ScrambledRadicalInverseSpecialized<5011>(perm, a);
    case 672:
        return ScrambledRadicalInverseSpecialized<5021>(perm, a);
    case 673:
        return ScrambledRadicalInverseSpecialized<5023>(perm, a);
    case 674:
        return ScrambledRadicalInverseSpecialized<5039>(perm, a);
    case 675:
        return ScrambledRadicalInverseSpecialized<5051>(perm, a);
    case 676:
        return ScrambledRadicalInverseSpecialized<5059>(perm, a);
    case 677:
        return ScrambledRadicalInverseSpecialized<5077>(perm, a);
    case 678:
        return ScrambledRadicalInverseSpecialized<5081>(perm, a);
    case 679:
        return ScrambledRadicalInverseSpecialized<5087>(perm, a);
    case 680:
        return ScrambledRadicalInverseSpecialized<5099>(perm, a);
    case 681:
        return ScrambledRadicalInverseSpecialized<5101>(perm, a);
    case 682:
        return ScrambledRadicalInverseSpecialized<5107>(perm, a);
    case 683:
        return ScrambledRadicalInverseSpecialized<5113>(perm, a);
    case 684:
        return ScrambledRadicalInverseSpecialized<5119>(perm, a);
    case 685:
        return ScrambledRadicalInverseSpecialized<5147>(perm, a);
    case 686:
        return ScrambledRadicalInverseSpecialized<5153>(perm, a);
    case 687:
        return ScrambledRadicalInverseSpecialized<5167>(perm, a);
    case 688:
        return ScrambledRadicalInverseSpecialized<5171>(perm, a);
    case 689:
        return ScrambledRadicalInverseSpecialized<5179>(perm, a);
    case 690:
        return ScrambledRadicalInverseSpecialized<5189>(perm, a);
    case 691:
        return ScrambledRadicalInverseSpecialized<5197>(perm, a);
    case 692:
        return ScrambledRadicalInverseSpecialized<5209>(perm, a);
    case 693:
        return ScrambledRadicalInverseSpecialized<5227>(perm, a);
    case 694:
        return ScrambledRadicalInverseSpecialized<5231>(perm, a);
    case 695:
        return ScrambledRadicalInverseSpecialized<5233>(perm, a);
    case 696:
        return ScrambledRadicalInverseSpecialized<5237>(perm, a);
    case 697:
        return ScrambledRadicalInverseSpecialized<5261>(perm, a);
    case 698:
        return ScrambledRadicalInverseSpecialized<5273>(perm, a);
    case 699:
        return ScrambledRadicalInverseSpecialized<5279>(perm, a);
    case 700:
        return ScrambledRadicalInverseSpecialized<5281>(perm, a);
    case 701:
        return ScrambledRadicalInverseSpecialized<5297>(perm, a);
    case 702:
        return ScrambledRadicalInverseSpecialized<5303>(perm, a);
    case 703:
        return ScrambledRadicalInverseSpecialized<5309>(perm, a);
    case 704:
        return ScrambledRadicalInverseSpecialized<5323>(perm, a);
    case 705:
        return ScrambledRadicalInverseSpecialized<5333>(perm, a);
    case 706:
        return ScrambledRadicalInverseSpecialized<5347>(perm, a);
    case 707:
        return ScrambledRadicalInverseSpecialized<5351>(perm, a);
    case 708:
        return ScrambledRadicalInverseSpecialized<5381>(perm, a);
    case 709:
        return ScrambledRadicalInverseSpecialized<5387>(perm, a);
    case 710:
        return ScrambledRadicalInverseSpecialized<5393>(perm, a);
    case 711:
        return ScrambledRadicalInverseSpecialized<5399>(perm, a);
    case 712:
        return ScrambledRadicalInverseSpecialized<5407>(perm, a);
    case 713:
        return ScrambledRadicalInverseSpecialized<5413>(perm, a);
    case 714:
        return ScrambledRadicalInverseSpecialized<5417>(perm, a);
    case 715:
        return ScrambledRadicalInverseSpecialized<5419>(perm, a);
    case 716:
        return ScrambledRadicalInverseSpecialized<5431>(perm, a);
    case 717:
        return ScrambledRadicalInverseSpecialized<5437>(perm, a);
    case 718:
        return ScrambledRadicalInverseSpecialized<5441>(perm, a);
    case 719:
        return ScrambledRadicalInverseSpecialized<5443>(perm, a);
    case 720:
        return ScrambledRadicalInverseSpecialized<5449>(perm, a);
    case 721:
        return ScrambledRadicalInverseSpecialized<5471>(perm, a);
    case 722:
        return ScrambledRadicalInverseSpecialized<5477>(perm, a);
    case 723:
        return ScrambledRadicalInverseSpecialized<5479>(perm, a);
    case 724:
        return ScrambledRadicalInverseSpecialized<5483>(perm, a);
    case 725:
        return ScrambledRadicalInverseSpecialized<5501>(perm, a);
    case 726:
        return ScrambledRadicalInverseSpecialized<5503>(perm, a);
    case 727:
        return ScrambledRadicalInverseSpecialized<5507>(perm, a);
    case 728:
        return ScrambledRadicalInverseSpecialized<5519>(perm, a);
    case 729:
        return ScrambledRadicalInverseSpecialized<5521>(perm, a);
    case 730:
        return ScrambledRadicalInverseSpecialized<5527>(perm, a);
    case 731:
        return ScrambledRadicalInverseSpecialized<5531>(perm, a);
    case 732:
        return ScrambledRadicalInverseSpecialized<5557>(perm, a);
    case 733:
        return ScrambledRadicalInverseSpecialized<5563>(perm, a);
    case 734:
        return ScrambledRadicalInverseSpecialized<5569>(perm, a);
    case 735:
        return ScrambledRadicalInverseSpecialized<5573>(perm, a);
    case 736:
        return ScrambledRadicalInverseSpecialized<5581>(perm, a);
    case 737:
        return ScrambledRadicalInverseSpecialized<5591>(perm, a);
    case 738:
        return ScrambledRadicalInverseSpecialized<5623>(perm, a);
    case 739:
        return ScrambledRadicalInverseSpecialized<5639>(perm, a);
    case 740:
        return ScrambledRadicalInverseSpecialized<5641>(perm, a);
    case 741:
        return ScrambledRadicalInverseSpecialized<5647>(perm, a);
    case 742:
        return ScrambledRadicalInverseSpecialized<5651>(perm, a);
    case 743:
        return ScrambledRadicalInverseSpecialized<5653>(perm, a);
    case 744:
        return ScrambledRadicalInverseSpecialized<5657>(perm, a);
    case 745:
        return ScrambledRadicalInverseSpecialized<5659>(perm, a);
    case 746:
        return ScrambledRadicalInverseSpecialized<5669>(perm, a);
    case 747:
        return ScrambledRadicalInverseSpecialized<5683>(perm, a);
    case 748:
        return ScrambledRadicalInverseSpecialized<5689>(perm, a);
    case 749:
        return ScrambledRadicalInverseSpecialized<5693>(perm, a);
    case 750:
        return ScrambledRadicalInverseSpecialized<5701>(perm, a);
    case 751:
        return ScrambledRadicalInverseSpecialized<5711>(perm, a);
    case 752:
        return ScrambledRadicalInverseSpecialized<5717>(perm, a);
    case 753:
        return ScrambledRadicalInverseSpecialized<5737>(perm, a);
    case 754:
        return ScrambledRadicalInverseSpecialized<5741>(perm, a);
    case 755:
        return ScrambledRadicalInverseSpecialized<5743>(perm, a);
    case 756:
        return ScrambledRadicalInverseSpecialized<5749>(perm, a);
    case 757:
        return ScrambledRadicalInverseSpecialized<5779>(perm, a);
    case 758:
        return ScrambledRadicalInverseSpecialized<5783>(perm, a);
    case 759:
        return ScrambledRadicalInverseSpecialized<5791>(perm, a);
    case 760:
        return ScrambledRadicalInverseSpecialized<5801>(perm, a);
    case 761:
        return ScrambledRadicalInverseSpecialized<5807>(perm, a);
    case 762:
        return ScrambledRadicalInverseSpecialized<5813>(perm, a);
    case 763:
        return ScrambledRadicalInverseSpecialized<5821>(perm, a);
    case 764:
        return ScrambledRadicalInverseSpecialized<5827>(perm, a);
    case 765:
        return ScrambledRadicalInverseSpecialized<5839>(perm, a);
    case 766:
        return ScrambledRadicalInverseSpecialized<5843>(perm, a);
    case 767:
        return ScrambledRadicalInverseSpecialized<5849>(perm, a);
    case 768:
        return ScrambledRadicalInverseSpecialized<5851>(perm, a);
    case 769:
        return ScrambledRadicalInverseSpecialized<5857>(perm, a);
    case 770:
        return ScrambledRadicalInverseSpecialized<5861>(perm, a);
    case 771:
        return ScrambledRadicalInverseSpecialized<5867>(perm, a);
    case 772:
        return ScrambledRadicalInverseSpecialized<5869>(perm, a);
    case 773:
        return ScrambledRadicalInverseSpecialized<5879>(perm, a);
    case 774:
        return ScrambledRadicalInverseSpecialized<5881>(perm, a);
    case 775:
        return ScrambledRadicalInverseSpecialized<5897>(perm, a);
    case 776:
        return ScrambledRadicalInverseSpecialized<5903>(perm, a);
    case 777:
        return ScrambledRadicalInverseSpecialized<5923>(perm, a);
    case 778:
        return ScrambledRadicalInverseSpecialized<5927>(perm, a);
    case 779:
        return ScrambledRadicalInverseSpecialized<5939>(perm, a);
    case 780:
        return ScrambledRadicalInverseSpecialized<5953>(perm, a);
    case 781:
        return ScrambledRadicalInverseSpecialized<5981>(perm, a);
    case 782:
        return ScrambledRadicalInverseSpecialized<5987>(perm, a);
    case 783:
        return ScrambledRadicalInverseSpecialized<6007>(perm, a);
    case 784:
        return ScrambledRadicalInverseSpecialized<6011>(perm, a);
    case 785:
        return ScrambledRadicalInverseSpecialized<6029>(perm, a);
    case 786:
        return ScrambledRadicalInverseSpecialized<6037>(perm, a);
    case 787:
        return ScrambledRadicalInverseSpecialized<6043>(perm, a);
    case 788:
        return ScrambledRadicalInverseSpecialized<6047>(perm, a);
    case 789:
        return ScrambledRadicalInverseSpecialized<6053>(perm, a);
    case 790:
        return ScrambledRadicalInverseSpecialized<6067>(perm, a);
    case 791:
        return ScrambledRadicalInverseSpecialized<6073>(perm, a);
    case 792:
        return ScrambledRadicalInverseSpecialized<6079>(perm, a);
    case 793:
        return ScrambledRadicalInverseSpecialized<6089>(perm, a);
    case 794:
        return ScrambledRadicalInverseSpecialized<6091>(perm, a);
    case 795:
        return ScrambledRadicalInverseSpecialized<6101>(perm, a);
    case 796:
        return ScrambledRadicalInverseSpecialized<6113>(perm, a);
    case 797:
        return ScrambledRadicalInverseSpecialized<6121>(perm, a);
    case 798:
        return ScrambledRadicalInverseSpecialized<6131>(perm, a);
    case 799:
        return ScrambledRadicalInverseSpecialized<6133>(perm, a);
    case 800:
        return ScrambledRadicalInverseSpecialized<6143>(perm, a);
    case 801:
        return ScrambledRadicalInverseSpecialized<6151>(perm, a);
    case 802:
        return ScrambledRadicalInverseSpecialized<6163>(perm, a);
    case 803:
        return ScrambledRadicalInverseSpecialized<6173>(perm, a);
    case 804:
        return ScrambledRadicalInverseSpecialized<6197>(perm, a);
    case 805:
        return ScrambledRadicalInverseSpecialized<6199>(perm, a);
    case 806:
        return ScrambledRadicalInverseSpecialized<6203>(perm, a);
    case 807:
        return ScrambledRadicalInverseSpecialized<6211>(perm, a);
    case 808:
        return ScrambledRadicalInverseSpecialized<6217>(perm, a);
    case 809:
        return ScrambledRadicalInverseSpecialized<6221>(perm, a);
    case 810:
        return ScrambledRadicalInverseSpecialized<6229>(perm, a);
    case 811:
        return ScrambledRadicalInverseSpecialized<6247>(perm, a);
    case 812:
        return ScrambledRadicalInverseSpecialized<6257>(perm, a);
    case 813:
        return ScrambledRadicalInverseSpecialized<6263>(perm, a);
    case 814:
        return ScrambledRadicalInverseSpecialized<6269>(perm, a);
    case 815:
        return ScrambledRadicalInverseSpecialized<6271>(perm, a);
    case 816:
        return ScrambledRadicalInverseSpecialized<6277>(perm, a);
    case 817:
        return ScrambledRadicalInverseSpecialized<6287>(perm, a);
    case 818:
        return ScrambledRadicalInverseSpecialized<6299>(perm, a);
    case 819:
        return ScrambledRadicalInverseSpecialized<6301>(perm, a);
    case 820:
        return ScrambledRadicalInverseSpecialized<6311>(perm, a);
    case 821:
        return ScrambledRadicalInverseSpecialized<6317>(perm, a);
    case 822:
        return ScrambledRadicalInverseSpecialized<6323>(perm, a);
    case 823:
        return ScrambledRadicalInverseSpecialized<6329>(perm, a);
    case 824:
        return ScrambledRadicalInverseSpecialized<6337>(perm, a);
    case 825:
        return ScrambledRadicalInverseSpecialized<6343>(perm, a);
    case 826:
        return ScrambledRadicalInverseSpecialized<6353>(perm, a);
    case 827:
        return ScrambledRadicalInverseSpecialized<6359>(perm, a);
    case 828:
        return ScrambledRadicalInverseSpecialized<6361>(perm, a);
    case 829:
        return ScrambledRadicalInverseSpecialized<6367>(perm, a);
    case 830:
        return ScrambledRadicalInverseSpecialized<6373>(perm, a);
    case 831:
        return ScrambledRadicalInverseSpecialized<6379>(perm, a);
    case 832:
        return ScrambledRadicalInverseSpecialized<6389>(perm, a);
    case 833:
        return ScrambledRadicalInverseSpecialized<6397>(perm, a);
    case 834:
        return ScrambledRadicalInverseSpecialized<6421>(perm, a);
    case 835:
        return ScrambledRadicalInverseSpecialized<6427>(perm, a);
    case 836:
        return ScrambledRadicalInverseSpecialized<6449>(perm, a);
    case 837:
        return ScrambledRadicalInverseSpecialized<6451>(perm, a);
    case 838:
        return ScrambledRadicalInverseSpecialized<6469>(perm, a);
    case 839:
        return ScrambledRadicalInverseSpecialized<6473>(perm, a);
    case 840:
        return ScrambledRadicalInverseSpecialized<6481>(perm, a);
    case 841:
        return ScrambledRadicalInverseSpecialized<6491>(perm, a);
    case 842:
        return ScrambledRadicalInverseSpecialized<6521>(perm, a);
    case 843:
        return ScrambledRadicalInverseSpecialized<6529>(perm, a);
    case 844:
        return ScrambledRadicalInverseSpecialized<6547>(perm, a);
    case 845:
        return ScrambledRadicalInverseSpecialized<6551>(perm, a);
    case 846:
        return ScrambledRadicalInverseSpecialized<6553>(perm, a);
    case 847:
        return ScrambledRadicalInverseSpecialized<6563>(perm, a);
    case 848:
        return ScrambledRadicalInverseSpecialized<6569>(perm, a);
    case 849:
        return ScrambledRadicalInverseSpecialized<6571>(perm, a);
    case 850:
        return ScrambledRadicalInverseSpecialized<6577>(perm, a);
    case 851:
        return ScrambledRadicalInverseSpecialized<6581>(perm, a);
    case 852:
        return ScrambledRadicalInverseSpecialized<6599>(perm, a);
    case 853:
        return ScrambledRadicalInverseSpecialized<6607>(perm, a);
    case 854:
        return ScrambledRadicalInverseSpecialized<6619>(perm, a);
    case 855:
        return ScrambledRadicalInverseSpecialized<6637>(perm, a);
    case 856:
        return ScrambledRadicalInverseSpecialized<6653>(perm, a);
    case 857:
        return ScrambledRadicalInverseSpecialized<6659>(perm, a);
    case 858:
        return ScrambledRadicalInverseSpecialized<6661>(perm, a);
    case 859:
        return ScrambledRadicalInverseSpecialized<6673>(perm, a);
    case 860:
        return ScrambledRadicalInverseSpecialized<6679>(perm, a);
    case 861:
        return ScrambledRadicalInverseSpecialized<6689>(perm, a);
    case 862:
        return ScrambledRadicalInverseSpecialized<6691>(perm, a);
    case 863:
        return ScrambledRadicalInverseSpecialized<6701>(perm, a);
    case 864:
        return ScrambledRadicalInverseSpecialized<6703>(perm, a);
    case 865:
        return ScrambledRadicalInverseSpecialized<6709>(perm, a);
    case 866:
        return ScrambledRadicalInverseSpecialized<6719>(perm, a);
    case 867:
        return ScrambledRadicalInverseSpecialized<6733>(perm, a);
    case 868:
        return ScrambledRadicalInverseSpecialized<6737>(perm, a);
    case 869:
        return ScrambledRadicalInverseSpecialized<6761>(perm, a);
    case 870:
        return ScrambledRadicalInverseSpecialized<6763>(perm, a);
    case 871:
        return ScrambledRadicalInverseSpecialized<6779>(perm, a);
    case 872:
        return ScrambledRadicalInverseSpecialized<6781>(perm, a);
    case 873:
        return ScrambledRadicalInverseSpecialized<6791>(perm, a);
    case 874:
        return ScrambledRadicalInverseSpecialized<6793>(perm, a);
    case 875:
        return ScrambledRadicalInverseSpecialized<6803>(perm, a);
    case 876:
        return ScrambledRadicalInverseSpecialized<6823>(perm, a);
    case 877:
        return ScrambledRadicalInverseSpecialized<6827>(perm, a);
    case 878:
        return ScrambledRadicalInverseSpecialized<6829>(perm, a);
    case 879:
        return ScrambledRadicalInverseSpecialized<6833>(perm, a);
    case 880:
        return ScrambledRadicalInverseSpecialized<6841>(perm, a);
    case 881:
        return ScrambledRadicalInverseSpecialized<6857>(perm, a);
    case 882:
        return ScrambledRadicalInverseSpecialized<6863>(perm, a);
    case 883:
        return ScrambledRadicalInverseSpecialized<6869>(perm, a);
    case 884:
        return ScrambledRadicalInverseSpecialized<6871>(perm, a);
    case 885:
        return ScrambledRadicalInverseSpecialized<6883>(perm, a);
    case 886:
        return ScrambledRadicalInverseSpecialized<6899>(perm, a);
    case 887:
        return ScrambledRadicalInverseSpecialized<6907>(perm, a);
    case 888:
        return ScrambledRadicalInverseSpecialized<6911>(perm, a);
    case 889:
        return ScrambledRadicalInverseSpecialized<6917>(perm, a);
    case 890:
        return ScrambledRadicalInverseSpecialized<6947>(perm, a);
    case 891:
        return ScrambledRadicalInverseSpecialized<6949>(perm, a);
    case 892:
        return ScrambledRadicalInverseSpecialized<6959>(perm, a);
    case 893:
        return ScrambledRadicalInverseSpecialized<6961>(perm, a);
    case 894:
        return ScrambledRadicalInverseSpecialized<6967>(perm, a);
    case 895:
        return ScrambledRadicalInverseSpecialized<6971>(perm, a);
    case 896:
        return ScrambledRadicalInverseSpecialized<6977>(perm, a);
    case 897:
        return ScrambledRadicalInverseSpecialized<6983>(perm, a);
    case 898:
        return ScrambledRadicalInverseSpecialized<6991>(perm, a);
    case 899:
        return ScrambledRadicalInverseSpecialized<6997>(perm, a);
    case 900:
        return ScrambledRadicalInverseSpecialized<7001>(perm, a);
    case 901:
        return ScrambledRadicalInverseSpecialized<7013>(perm, a);
    case 902:
        return ScrambledRadicalInverseSpecialized<7019>(perm, a);
    case 903:
        return ScrambledRadicalInverseSpecialized<7027>(perm, a);
    case 904:
        return ScrambledRadicalInverseSpecialized<7039>(perm, a);
    case 905:
        return ScrambledRadicalInverseSpecialized<7043>(perm, a);
    case 906:
        return ScrambledRadicalInverseSpecialized<7057>(perm, a);
    case 907:
        return ScrambledRadicalInverseSpecialized<7069>(perm, a);
    case 908:
        return ScrambledRadicalInverseSpecialized<7079>(perm, a);
    case 909:
        return ScrambledRadicalInverseSpecialized<7103>(perm, a);
    case 910:
        return ScrambledRadicalInverseSpecialized<7109>(perm, a);
    case 911:
        return ScrambledRadicalInverseSpecialized<7121>(perm, a);
    case 912:
        return ScrambledRadicalInverseSpecialized<7127>(perm, a);
    case 913:
        return ScrambledRadicalInverseSpecialized<7129>(perm, a);
    case 914:
        return ScrambledRadicalInverseSpecialized<7151>(perm, a);
    case 915:
        return ScrambledRadicalInverseSpecialized<7159>(perm, a);
    case 916:
        return ScrambledRadicalInverseSpecialized<7177>(perm, a);
    case 917:
        return ScrambledRadicalInverseSpecialized<7187>(perm, a);
    case 918:
        return ScrambledRadicalInverseSpecialized<7193>(perm, a);
    case 919:
        return ScrambledRadicalInverseSpecialized<7207>(perm, a);
    case 920:
        return ScrambledRadicalInverseSpecialized<7211>(perm, a);
    case 921:
        return ScrambledRadicalInverseSpecialized<7213>(perm, a);
    case 922:
        return ScrambledRadicalInverseSpecialized<7219>(perm, a);
    case 923:
        return ScrambledRadicalInverseSpecialized<7229>(perm, a);
    case 924:
        return ScrambledRadicalInverseSpecialized<7237>(perm, a);
    case 925:
        return ScrambledRadicalInverseSpecialized<7243>(perm, a);
    case 926:
        return ScrambledRadicalInverseSpecialized<7247>(perm, a);
    case 927:
        return ScrambledRadicalInverseSpecialized<7253>(perm, a);
    case 928:
        return ScrambledRadicalInverseSpecialized<7283>(perm, a);
    case 929:
        return ScrambledRadicalInverseSpecialized<7297>(perm, a);
    case 930:
        return ScrambledRadicalInverseSpecialized<7307>(perm, a);
    case 931:
        return ScrambledRadicalInverseSpecialized<7309>(perm, a);
    case 932:
        return ScrambledRadicalInverseSpecialized<7321>(perm, a);
    case 933:
        return ScrambledRadicalInverseSpecialized<7331>(perm, a);
    case 934:
        return ScrambledRadicalInverseSpecialized<7333>(perm, a);
    case 935:
        return ScrambledRadicalInverseSpecialized<7349>(perm, a);
    case 936:
        return ScrambledRadicalInverseSpecialized<7351>(perm, a);
    case 937:
        return ScrambledRadicalInverseSpecialized<7369>(perm, a);
    case 938:
        return ScrambledRadicalInverseSpecialized<7393>(perm, a);
    case 939:
        return ScrambledRadicalInverseSpecialized<7411>(perm, a);
    case 940:
        return ScrambledRadicalInverseSpecialized<7417>(perm, a);
    case 941:
        return ScrambledRadicalInverseSpecialized<7433>(perm, a);
    case 942:
        return ScrambledRadicalInverseSpecialized<7451>(perm, a);
    case 943:
        return ScrambledRadicalInverseSpecialized<7457>(perm, a);
    case 944:
        return ScrambledRadicalInverseSpecialized<7459>(perm, a);
    case 945:
        return ScrambledRadicalInverseSpecialized<7477>(perm, a);
    case 946:
        return ScrambledRadicalInverseSpecialized<7481>(perm, a);
    case 947:
        return ScrambledRadicalInverseSpecialized<7487>(perm, a);
    case 948:
        return ScrambledRadicalInverseSpecialized<7489>(perm, a);
    case 949:
        return ScrambledRadicalInverseSpecialized<7499>(perm, a);
    case 950:
        return ScrambledRadicalInverseSpecialized<7507>(perm, a);
    case 951:
        return ScrambledRadicalInverseSpecialized<7517>(perm, a);
    case 952:
        return ScrambledRadicalInverseSpecialized<7523>(perm, a);
    case 953:
        return ScrambledRadicalInverseSpecialized<7529>(perm, a);
    case 954:
        return ScrambledRadicalInverseSpecialized<7537>(perm, a);
    case 955:
        return ScrambledRadicalInverseSpecialized<7541>(perm, a);
    case 956:
        return ScrambledRadicalInverseSpecialized<7547>(perm, a);
    case 957:
        return ScrambledRadicalInverseSpecialized<7549>(perm, a);
    case 958:
        return ScrambledRadicalInverseSpecialized<7559>(perm, a);
    case 959:
        return ScrambledRadicalInverseSpecialized<7561>(perm, a);
    case 960:
        return ScrambledRadicalInverseSpecialized<7573>(perm, a);
    case 961:
        return ScrambledRadicalInverseSpecialized<7577>(perm, a);
    case 962:
        return ScrambledRadicalInverseSpecialized<7583>(perm, a);
    case 963:
        return ScrambledRadicalInverseSpecialized<7589>(perm, a);
    case 964:
        return ScrambledRadicalInverseSpecialized<7591>(perm, a);
    case 965:
        return ScrambledRadicalInverseSpecialized<7603>(perm, a);
    case 966:
        return ScrambledRadicalInverseSpecialized<7607>(perm, a);
    case 967:
        return ScrambledRadicalInverseSpecialized<7621>(perm, a);
    case 968:
        return ScrambledRadicalInverseSpecialized<7639>(perm, a);
    case 969:
        return ScrambledRadicalInverseSpecialized<7643>(perm, a);
    case 970:
        return ScrambledRadicalInverseSpecialized<7649>(perm, a);
    case 971:
        return ScrambledRadicalInverseSpecialized<7669>(perm, a);
    case 972:
        return ScrambledRadicalInverseSpecialized<7673>(perm, a);
    case 973:
        return ScrambledRadicalInverseSpecialized<7681>(perm, a);
    case 974:
        return ScrambledRadicalInverseSpecialized<7687>(perm, a);
    case 975:
        return ScrambledRadicalInverseSpecialized<7691>(perm, a);
    case 976:
        return ScrambledRadicalInverseSpecialized<7699>(perm, a);
    case 977:
        return ScrambledRadicalInverseSpecialized<7703>(perm, a);
    case 978:
        return ScrambledRadicalInverseSpecialized<7717>(perm, a);
    case 979:
        return ScrambledRadicalInverseSpecialized<7723>(perm, a);
    case 980:
        return ScrambledRadicalInverseSpecialized<7727>(perm, a);
    case 981:
        return ScrambledRadicalInverseSpecialized<7741>(perm, a);
    case 982:
        return ScrambledRadicalInverseSpecialized<7753>(perm, a);
    case 983:
        return ScrambledRadicalInverseSpecialized<7757>(perm, a);
    case 984:
        return ScrambledRadicalInverseSpecialized<7759>(perm, a);
    case 985:
        return ScrambledRadicalInverseSpecialized<7789>(perm, a);
    case 986:
        return ScrambledRadicalInverseSpecialized<7793>(perm, a);
    case 987:
        return ScrambledRadicalInverseSpecialized<7817>(perm, a);
    case 988:
        return ScrambledRadicalInverseSpecialized<7823>(perm, a);
    case 989:
        return ScrambledRadicalInverseSpecialized<7829>(perm, a);
    case 990:
        return ScrambledRadicalInverseSpecialized<7841>(perm, a);
    case 991:
        return ScrambledRadicalInverseSpecialized<7853>(perm, a);
    case 992:
        return ScrambledRadicalInverseSpecialized<7867>(perm, a);
    case 993:
        return ScrambledRadicalInverseSpecialized<7873>(perm, a);
    case 994:
        return ScrambledRadicalInverseSpecialized<7877>(perm, a);
    case 995:
        return ScrambledRadicalInverseSpecialized<7879>(perm, a);
    case 996:
        return ScrambledRadicalInverseSpecialized<7883>(perm, a);
    case 997:
        return ScrambledRadicalInverseSpecialized<7901>(perm, a);
    case 998:
        return ScrambledRadicalInverseSpecialized<7907>(perm, a);
    case 999:
        return ScrambledRadicalInverseSpecialized<7919>(perm, a);
    case 1000:
        return ScrambledRadicalInverseSpecialized<7927>(perm, a);
    case 1001:
        return ScrambledRadicalInverseSpecialized<7933>(perm, a);
    case 1002:
        return ScrambledRadicalInverseSpecialized<7937>(perm, a);
    case 1003:
        return ScrambledRadicalInverseSpecialized<7949>(perm, a);
    case 1004:
        return ScrambledRadicalInverseSpecialized<7951>(perm, a);
    case 1005:
        return ScrambledRadicalInverseSpecialized<7963>(perm, a);
    case 1006:
        return ScrambledRadicalInverseSpecialized<7993>(perm, a);
    case 1007:
        return ScrambledRadicalInverseSpecialized<8009>(perm, a);
    case 1008:
        return ScrambledRadicalInverseSpecialized<8011>(perm, a);
    case 1009:
        return ScrambledRadicalInverseSpecialized<8017>(perm, a);
    case 1010:
        return ScrambledRadicalInverseSpecialized<8039>(perm, a);
    case 1011:
        return ScrambledRadicalInverseSpecialized<8053>(perm, a);
    case 1012:
        return ScrambledRadicalInverseSpecialized<8059>(perm, a);
    case 1013:
        return ScrambledRadicalInverseSpecialized<8069>(perm, a);
    case 1014:
        return ScrambledRadicalInverseSpecialized<8081>(perm, a);
    case 1015:
        return ScrambledRadicalInverseSpecialized<8087>(perm, a);
    case 1016:
        return ScrambledRadicalInverseSpecialized<8089>(perm, a);
    case 1017:
        return ScrambledRadicalInverseSpecialized<8093>(perm, a);
    case 1018:
        return ScrambledRadicalInverseSpecialized<8101>(perm, a);
    case 1019:
        return ScrambledRadicalInverseSpecialized<8111>(perm, a);
    case 1020:
        return ScrambledRadicalInverseSpecialized<8117>(perm, a);
    case 1021:
        return ScrambledRadicalInverseSpecialized<8123>(perm, a);
    case 1022:
        return ScrambledRadicalInverseSpecialized<8147>(perm, a);
    case 1023:
        return ScrambledRadicalInverseSpecialized<8161>(perm, a);
    default:
        // LOG(FATAL) << StringPrintf("Base %d is >= 1024, the limit of "
        //                            "ScrambledRadicalInverse",
        //                            baseIndex);
        return 0;
    }
}


} // namespace pbrt