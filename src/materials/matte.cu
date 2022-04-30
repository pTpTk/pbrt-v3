
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


// materials/matte.cpp*
#include "materials/matte.cuh"
#include "paramset.cuh"
#include "reflection.cuh"
#include "interaction.cuh"
#include "texture.cuh"
#include "interaction.cuh"

namespace pbrt {

// MatteMaterial Method Definitions
__both__
void MatteMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                               MemoryArena &arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const {
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);

    // Evaluate textures for _MatteMaterial_ material and allocate BRDF
    si->bsdf = new BSDF(*si);
    Spectrum r = Kd->Evaluate(*si).Clamp();
    Float sig = Clamp(sigma->Evaluate(*si), 0, 90);
    if (!r.IsBlack()) {
        if (sig == 0)
            si->bsdf->Add(new LambertianReflection(r));
    }
}

MatteMaterial *CreateMatteMaterial(const TextureParams &mp) {
    Texture<Spectrum>* Kd;// = mp.GetSpectrumTexture("Kd", Spectrum(0.5f)).get();
    Texture<Float>* sigma;// = mp.GetFloatTexture("sigma", 0.f).get();
    Texture<Float>* bumpMap;// = mp.GetFloatTextureOrNull("bumpmap").get();
    
    cudaMallocManaged(&Kd, sizeof(Texture<Spectrum>));
    std::cout << "Error matte.cu: 69: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    cudaMallocManaged(&sigma, sizeof(Texture<Float>));
    std::cout << "Error matte.cu: 71: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    cudaMallocManaged(&bumpMap, sizeof(Texture<Float>));
    std::cout << "Error matte.cu: 73: " << cudaGetErrorString(cudaGetLastError()) << std::endl;

    *Kd = *mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    *sigma = *mp.GetFloatTexture("sigma", 0.f);
    *bumpMap = *mp.GetFloatTextureOrNull("bumpmap");

    void* ptr;
    cudaMallocManaged(&ptr, sizeof(MatteMaterial));
    std::cout << "Error matte.cu: 81: " << cudaGetErrorString(cudaGetLastError()) << std::endl;
    return new(ptr) MatteMaterial(Kd, sigma, bumpMap);
}

}  // namespace pbrt
