// accelerators/accel-bvh.cpp*
#include "accelerators/accel-bvh.h"

namespace pbrt {

AccelBVHAccel::AccelBVHAccel(std::vector<std::shared_ptr<Primitive>> p) {
    std::cout << "This is AccelBVHAccel::AccelBVHAccel\n";
    assert(0 == -1);
}

Bounds3f AccelBVHAccel::WorldBound() const {}

bool AccelBVHAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {}
bool AccelBVHAccel::IntersectP(const Ray &ray) const {}


std::shared_ptr<AccelBVHAccel> CreateAccelBVHAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims) {
    return std::make_shared<AccelBVHAccel>(std::move(prims));
}

} // namespace pbrt
