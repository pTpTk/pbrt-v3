#include "gpu.cuh"

namespace pbrt{
namespace gpu{
    
__device__
Spectrum Li(const RayDifferential &r, const Scene &scene, Sampler &sampler, int depth){
    Spectrum L(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
    int bounces;
    Float etaScale = 1;

    for (bounces = 0;; ++bounces) {
    
        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

    }
}

};
};
