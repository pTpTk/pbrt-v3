#ifndef PBRT_ACCELERATORS_ACCEL_BVH_H
#define PBRT_ACCELERATORS_ACCEL_BVH_H

// accelerators/accel-bvh.h*
#include "pbrt.h"
#include "primitive.h"
#include <atomic>

namespace pbrt {
class AccelBVHAccel : public Aggregate {
  public:
    AccelBVHAccel(std::vector<std::shared_ptr<Primitive>> p);

    Bounds3f WorldBound() const;
    bool Intersect(const Ray &r, SurfaceInteraction *) const;
    bool IntersectP(const Ray &ray) const;
  private:

};

std::shared_ptr<AccelBVHAccel> CreateAccelBVHAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims);

} // namespace pbrt

#endif // PBRT_ACCELERATORS_ACCEL_BVH_H
