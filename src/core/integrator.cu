
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

// core/integrator.cpp*
#include "integrator.cuh"
#include "scene.cuh"
#include "interaction.cuh"
#include "sampling.cuh"
#include "parallel.cuh"
#include "film.cuh"
#include "sampler.cuh"
#include "integrator.cuh"
#include "progressreporter.cuh"
#include "camera.cuh"
#include "stats.cuh"

namespace pbrt {



// Integrator Method Definitions
Integrator::~Integrator() {
    
}

// Integrator Utility Functions



// SamplerIntegrator Method Definitions
// void SamplerIntegrator::Render(Scene *s) {
//     Scene& scene = *s;
//     CallLiKernel();
//     Preprocess(scene, *sampler);
//     // Render image tiles in parallel

//     // Compute number of tiles, _nTiles_, to use for parallel rendering
//     Bounds2i sampleBounds = camera->film->GetSampleBounds();
//     Vector2i sampleExtent = sampleBounds.Diagonal();
//     const int tileSize = 16;
//     Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
//                    (sampleExtent.y + tileSize - 1) / tileSize);
//     LOG(ERROR) << "\n" << cudaGetErrorString(cudaGetLastError()) << std::endl;
//     CallLiKernel();
//     // ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
//     // {
//     //     printf("here\n");
//     //     LiKernel<<<1,1>>>();
//     //     std::cout << "\n" << cudaGetErrorString(cudaGetLastError()) << std::endl;
//     //     // ParallelFor2D([&](Point2i tile) {
//     //     //     // Render section of image corresponding to _tile_

//     //     //     // Allocate _MemoryArena_ for tile
//     //     //     MemoryArena arena;

//     //     //     // Get sampler instance for tile
//     //     //     int seed = tile.y * nTiles.x + tile.x;
//     //     //     std::vector<Sampler*> tileSamplers;

//     //     //     // Compute sample bounds for tile
//     //     //     int x0 = sampleBounds.pMin.x + tile.x * tileSize;
//     //     //     int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
//     //     //     int y0 = sampleBounds.pMin.y + tile.y * tileSize;
//     //     //     int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
//     //     //     Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
//     //     //     LOG(INFO) << "Starting image tile " << tileBounds;

//     //     //     // Get _FilmTile_ for tile
//     //     //     std::unique_ptr<FilmTile> filmTile =
//     //     //         camera->film->GetFilmTile(tileBounds);

//     //     //     // Creating unified objects for calculations
//     //     //     CameraSample* cameraSamples;
//     //     //     RayDifferential* rays;
//     //     //     Float* rayWeights;
//     //     //     Spectrum* Ls;
//     //     //     cudaMallocManaged(&cameraSamples, sizeof(CameraSample) 
//     //     //                       * sampler->samplesPerPixel);
//     //     //     cudaMallocManaged(&rays,          sizeof(RayDifferential) 
//     //     //                       * sampler->samplesPerPixel); 
//     //     //     cudaMallocManaged(&rayWeights,    sizeof(Float) 
//     //     //                       * sampler->samplesPerPixel); 
//     //     //     cudaMallocManaged(&Ls,            sizeof(Spectrum) 
//     //     //                       * sampler->samplesPerPixel); 

//     //     //     new(cameraSamples) CameraSample[sampler->samplesPerPixel];
//     //     //     new(rays)         RayDifferential[sampler->samplesPerPixel];
//     //     //     new(rayWeights)   Float[sampler->samplesPerPixel];
//     //     //     new(Ls)           Spectrum[sampler->samplesPerPixel];

//     //     //     // Loop over pixels in tile to render them
//     //     //     for (Point2i pixel : tileBounds) {
//     //     //         if (!InsideExclusive(pixel, pixelBounds))
//     //     //             continue;

//     //     //         for(int64_t sampleNum = 0; 
//     //     //             sampleNum < sampler->samplesPerPixel; 
//     //     //             sampleNum++) {

//     //     //             tileSamplers.push_back(sampler->Clone(seed));
//     //     //             {
//     //     //                 ProfilePhase pp(Prof::StartPixel);
//     //     //                 tileSamplers[sampleNum]->StartPixel(pixel);
//     //     //             }
//     //     //             tileSamplers[sampleNum]->SetSampleNumber(sampleNum);
                
//     //     //             // Initialize _CameraSample_ for current sample
//     //     //             cameraSamples[sampleNum] =
//     //     //                 tileSamplers[sampleNum]->GetCameraSample(pixel);

//     //     //             // Generate camera ray for current sample
                    
//     //     //             rayWeights[sampleNum] =
//     //     //                 camera->GenerateRayDifferential(cameraSamples[sampleNum], &rays[sampleNum]);
//     //     //             rays[sampleNum].ScaleDifferentials(
//     //     //                 1 / std::sqrt((Float)tileSamplers[sampleNum]->samplesPerPixel));
//     //     //             ++nCameraRays;
//     //     //             // printf("Host: Ls[%d] address: %x\n", sampleNum, &Ls[sampleNum]);
//     //     //         }
                
//     //     //         // Evaluate radiance along camera ray
//     //     //         LiKernel<<<1, sampler->samplesPerPixel>>>
//     //     //             (Ls, this, rays, rayWeights, scene, tileSamplers.data(), arena);
//     //     //         cudaDeviceSynchronize();

//     //     //         for(int64_t sampleNum = 0; 
//     //     //             sampleNum < sampler->samplesPerPixel; 
//     //     //             sampleNum++) {
                    
//     //     //             // Ls[sampleNum] = 0.f;
//     //     //             // if (rayWeights[sampleNum] > 0) 
//     //     //             //     Ls[sampleNum] = 
//     //     //             //         Li(rays[sampleNum], scene, *tileSamplers[sampleNum], arena);

//     //     //             // Issue warning if unexpected radiance value returned
//     //     //             // if(Ls[sampleNum] != Spectrum(20.0)){
//     //     //             //     LOG(ERROR) << StringPrintf(
//     //     //             //         "Error for pixel (%d, %d), sample %d.",
//     //     //             //         pixel.x, pixel.y,
//     //     //             //         (int)tileSamplers[sampleNum]->CurrentSampleNumber())
//     //     //             //         <<" Ls[sampleNum] = " << Ls[sampleNum];
//     //     //             // }

//     //     //             if (Ls[sampleNum].HasNaNs()) {
//     //     //                 LOG(ERROR) << StringPrintf(
//     //     //                     "Not-a-number radiance value returned "
//     //     //                     "for pixel (%d, %d), sample %d. Setting to black.",
//     //     //                     pixel.x, pixel.y,
//     //     //                     (int)tileSamplers[sampleNum]->CurrentSampleNumber());
//     //     //                 Ls[sampleNum] = Spectrum(0.f);
//     //     //             } else if (Ls[sampleNum].y() < -1e-5) {
//     //     //                 LOG(ERROR) << StringPrintf(
//     //     //                     "Negative luminance value, %f, returned "
//     //     //                     "for pixel (%d, %d), sample %d. Setting to black.",
//     //     //                     Ls[sampleNum].y(), pixel.x, pixel.y,
//     //     //                     (int)tileSamplers[sampleNum]->CurrentSampleNumber());
//     //     //                 Ls[sampleNum] = Spectrum(0.f);
//     //     //             } else if (std::isinf(Ls[sampleNum].y())) {
//     //     //                   LOG(ERROR) << StringPrintf(
//     //     //                     "Infinite luminance value returned "
//     //     //                     "for pixel (%d, %d), sample %d. Setting to black.",
//     //     //                     pixel.x, pixel.y,
//     //     //                     (int)tileSamplers[sampleNum]->CurrentSampleNumber());
//     //     //                 Ls[sampleNum] = Spectrum(0.f);
//     //     //             }
//     //     //             VLOG(1) << "Camera sample: " << cameraSamples[sampleNum] 
//     //     //                     << " -> ray: " << rays[sampleNum]
//     //     //                     << " -> Ls[sampleNum] = " << Ls[sampleNum];

//     //     //             // Add camera ray's contribution to image
//     //     //             filmTile->AddSample(cameraSamples[sampleNum].pFilm, 
//     //     //                                 Ls[sampleNum], rayWeights[sampleNum]);

//     //     //             // Free _MemoryArena_ memory from computing image sample
//     //     //             // value
//     //     //             arena.Reset();
//     //     //         }
//     //     //     }
//     //     //     LOG(INFO) << "Finished image tile " << tileBounds;

//     //     //     // Merge image tile into _Film_
//     //     //     camera->film->MergeFilmTile(std::move(filmTile));
//     //     //     reporter.Update();
//     //     //     cudaFree(cameraSamples);
//     //     //     cudaFree(rays);
//     //     //     cudaFree(rayWeights);
//     //     //     cudaFree(Ls);
//     //     // }, nTiles);
//     //     reporter.Done();
//     // }
//     // LOG(INFO) << "Rendering finished";

//     // Save final image after rendering
//     camera->film->WriteImage();

// }

// __global__
// void LiKernel(Spectrum* Ls, SamplerIntegrator* integrator,
//               const RayDifferential* rays, const Float* rayWeights,
//               const Scene &scene, Sampler** tileSamplers,
//               MemoryArena &arena) {
//     int sampleNum = threadIdx.x;
//     Ls[sampleNum] = Spectrum(20.0); // 0.f;
//     printf("Device: Ls[%d] address: %x\n", sampleNum, &Ls[sampleNum]);
//     // if (rayWeights[sampleNum] > 0)
//     //     Ls[sampleNum] = integrator->Li(rays[sampleNum], scene, *tileSamplers[sampleNum], arena);
// }

// __global__
// void LiKernel(int i, int* j){
//     // printf("Hello World\n");
//     *j = i;
// }

// void Render(Integrator *i, Scene *s){




}  // namespace pbrt
