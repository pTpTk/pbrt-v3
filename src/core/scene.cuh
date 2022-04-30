
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

// core/scene.h*
#include "pbrt.cuh"
#include "geometry.cuh"
#include "primitive.cuh"
#include "light.cuh"

namespace pbrt {

// Scene Declarations
class Scene {
  public:
    // Scene Public Methods
    Scene(Primitive* aggregate,
          std::vector<Light*> lights_in)
        : lights_v(lights_in), aggregate(aggregate) {
        // Scene Constructor Implementation
        worldBound = aggregate->WorldBound();
        std::size_t count = 0;
        for (int i = 0; i < lights_v.size(); ++i) {
            Light* light = lights_v[i];
            printf("LIGHT*: %x\n", (uint64_t)light);
            light->Preprocess(*this);
            if (light->flags & (int)LightFlags::Infinite) {
                ++count;
            }
        }
        infiniteLights_size = count;
        // allocate memory
        void* ptr;
        cudaMallocManaged(&ptr, sizeof(Light*) * count);
        LOG(ERROR) << "\n" << cudaGetErrorString(cudaGetLastError()) << std::endl;
        infiniteLights = new(ptr) Light*[count];
        count = 0;
        for (int i = 0; i < lights_v.size(); i++) {
            Light* light = lights_v[i];
            // uncomment this line if you see runtime probs
            // light->Preprocess(*this);
            if (light->flags & (int)LightFlags::Infinite) {
                infiniteLights[count++] = light;
            }
        }

        lights = (Light**)lights_v.data();
        lights_size = lights_v.size();


    }
    __both__
    const Bounds3f &WorldBound() const { return worldBound; }
    __both__
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    __both__
    bool IntersectP(const Ray &ray) const;
    __both__
    bool IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
                     Spectrum *transmittance) const;

    // Scene Public Data
    std::vector<Light*> lights_v;
    Light** lights;
    // Store infinite light sources separately for cases where we only want
    // to loop over them.
    Light** infiniteLights;
    int lights_size, infiniteLights_size;

  private:
    // Scene Private Data
    Primitive* aggregate;
    Bounds3f worldBound;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SCENE_H
