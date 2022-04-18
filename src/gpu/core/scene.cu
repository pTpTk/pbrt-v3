#include "gpu.cuh"

namespace pbrt {
namespace gpu{

__device__
Scene::Scene(){}

__device__
const Ray &ray, SurfaceInteraction *isect) const {
    assert(ray.d != Vector3f(0,0,0));
    return aggregate->Intersect(ray, isect);
}

};
};