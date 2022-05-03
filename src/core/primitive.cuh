
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

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
#include <atomic>
#include "accelerators/bvh.cuh"
#include "pbrt.cuh"
#include "shape.cuh"
#include "material.cuh"
#include "medium.cuh"
#include "transform.cuh"

#include <cstdint>

namespace pbrt {

struct BVHBuildNode;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct MortonPrimitive;
struct LinearBVHNode;

Primitive* CreateBVHAccelerator(
    std::vector<Primitive*> prims, const ParamSet &ps);

// Primitive Declarations
class Primitive {
  public:
    enum class PrimitiveType: std::uint8_t {
      Aggregate,
      GeometricPrimitive,
      BVHAccel,
      Primitive,
      Unspecified, // should never be this type
    };
  
    // Primitive Interface
    ~Primitive();
    __both__
    Bounds3f WorldBound() const;
    __both__
    bool Intersect(const Ray &r, SurfaceInteraction *) const;
    __both__
    bool IntersectP(const Ray &r) const;
    __both__
    const Light *GetAreaLight() const;
    __both__
    const Material *GetMaterial() const;
    __both__
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

    // GeometricPrimitive Interface
    // Set type to PrimitiveType::GeometricPrimitive in this constructor
    __both__
    Primitive(Shape*     const shape,
              Material*  const material,
              Light* const areaLight,
              const MediumInterface &mediumInterface);

    // BVHAccel Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    // BVHAccel Public Methods
    // Set type to PrimitiveType::BVHAccel in this constructor
    Primitive(std::vector<Primitive*> p,
              int maxPrimsInNode = 1,
              SplitMethod splitMethod = SplitMethod::SAH);

  private:
    // Primitive Private Data
    PrimitiveType type = PrimitiveType::Primitive;

    // GeometricPrimitive Private Data
    Shape* shape;
    Material* material;
    Light* areaLight;
    MediumInterface mediumInterface;

    // BVHAccel Private Methods
    BVHBuildNode *recursiveBuild(
        MemoryArena &arena, std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int start, int end, int *totalNodes,
        std::vector<Primitive*> &orderedPrims);
    BVHBuildNode *HLBVHBuild(
        MemoryArena &arena, const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        int *totalNodes,
        std::vector<Primitive*> &orderedPrims) const;
    BVHBuildNode *emitLBVH(
        BVHBuildNode *&buildNodes,
        const std::vector<BVHPrimitiveInfo> &primitiveInfo,
        MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
        std::vector<Primitive*> &orderedPrims,
        std::atomic<int> *orderedPrimsOffset, int bitIndex) const;
    BVHBuildNode *buildUpperSAH(MemoryArena &arena,
                                std::vector<BVHBuildNode *> &treeletRoots,
                                int start, int end, int *totalNodes) const;
    int flattenBVHTree(BVHBuildNode *node, int *offset);

    // BVHAccel Private Data
    int maxPrimsInNode;
    SplitMethod splitMethod;
    Primitive** primitives;
    LinearBVHNode *nodes = nullptr;
    std::vector<Primitive*> primitives_v;
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
