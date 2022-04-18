#include "gpu.cuh"
#include "geometry.h"
#include "interaction.h"
#include "shapes/sphere.cu"

namespace pbrt {
namespace gpu {

__device__
Shape::Shape(void* p, ShapeType type)
: p.p(p),
  type(type)
{}

__device__
bool Shape::Intersect(const Ray &ray, Float *tHit,
                      SurfaceInteraction *isect,
                      bool testAlphaTexture = true) const {
    switch(type){
        case ShapeType::SphereShape:
            return p.sphere->Intersect(ray, tHit, isect, testAlphaTexture);
        default:
            printf("Unknown Shape\n");
            assert(false);
    }
}

}
}