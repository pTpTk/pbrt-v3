#include "gpu.cuh"

namespace pbrt {
namespace gpu{

class Scene {
  public:
    __device__ Scene();
    __device__ bool Intersect(const Ray, SurfaceInteraction*) const;

  private:
    Primitive* aggregate;
};

};
};