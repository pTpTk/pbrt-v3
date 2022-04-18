#include "gpu.cuh"
#include "geometry.h"
#include "interaction.h"
#include "shapes/sphere.cu"

namespace pbrt {
namespace gpu {

union ShapePtr{
    void* p;
    Sphere* sphere;
};

enum struct ShapeType { SphereShape, TriangleShape };

class Shape {
  public:
    __device__ Shape(void*, ShapeType);
    __device__ bool Intersect(const Ray &ray, Float *tHit,
                              SurfaceInteraction *isect,
                              bool testAlphaTexture) const;
  private:
    ShapePtr p;
    ShapeType type;
};

}
}