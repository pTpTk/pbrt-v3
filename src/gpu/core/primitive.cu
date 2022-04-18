#include "gpu.cuh"
#include "shape.cu"
#include "geometry.h"
#include "interaction.cu"
#include "medium.cu"

namespace pbrt {
namespace gpu {

__device__
GeometricPrimitive::GeometricPrimitive(Shape* shape)
: shape(shape)
{}

__device__
bool GeometricPrimitive::Intersect(const Ray &r,
                              SurfaceInteraction *isect) const {
    Float tHit;
    if (!shape->Intersect(r, &tHit, isect)) return false;
    r.tMax = tHit;
    isect->primitive = this;
    assert(Dot(isect->n, isect->shading.n) => 0.);
    // Initialize _SurfaceInteraction::mediumInterface_ after _Shape_
    // intersection
    if (mediumInterface.IsMediumTransition())
        isect->mediumInterface = mediumInterface;
    else
        isect->mediumInterface = MediumInterface(r.medium);
    return true;

}


};
};