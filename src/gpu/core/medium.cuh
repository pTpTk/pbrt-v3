#include "gpu.cuh"

namespace pbrt {
namespace gpu {

class Medium {
  public:
    // Medium Interface
    __device__ virtual ~Medium() {}
    __device__ virtual Spectrum Tr(const Ray &ray, Sampler &sampler) const = 0;
    __device__ virtual Spectrum Sample(const Ray &ray, Sampler &sampler,
                            MemoryArena &arena,
                            MediumInteraction *mi) const = 0;
};

struct MediumInterface {
    __device__
    MediumInterface() : inside(nullptr), outside(nullptr) {}
    // MediumInterface Public Methods
    __device__
    MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}
    __device__
    MediumInterface(const Medium *inside, const Medium *outside)
        : inside(inside), outside(outside) {}
    __device__
    bool IsMediumTransition() const { return inside != outside; }
    const Medium *inside, *outside;
};

};
};