#include "gpu.cuh"
#include "shape.cu"
#include "geometry.h"
#include "interaction.cu"
#include "medium.cu"

namespace pbrt {
namespace gpu {

class GeometricPrimitive{
  public:
    __device__ GeometricPrimitive(Shape*);
    __device__ bool Intersect(const Ray &r, SurfaceInteraction *isect) const;

  private:
    Shape* shape;
};

};
};