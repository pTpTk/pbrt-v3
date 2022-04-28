
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

// core/api.cpp*
#include "api.cuh"
#include "parallel.cuh"
#include "paramset.cuh"
#include "spectrum.cuh"
#include "scene.cuh"
#include "film.cuh"
#include "medium.cuh"
#include "stats.cuh"

// API Additional Headers
#include "accelerators/bvh.cuh"
#include "cameras/perspective.cuh"
#include "filters/box.cuh"
#include "integrators/path.cuh"
#include "lights/diffuse.cuh"
#include "materials/matte.cuh"
#include "samplers/halton.cuh"
#include "shapes/sphere.cuh"

#include <map>
#include <stdio.h>

namespace pbrt {

// API Global Variables
Options PbrtOptions;

// API Local Classes
PBRT_CONSTEXPR int MaxTransforms = 2;
// PBRT_CONSTEXPR int StartTransformBits = 1 << 0;
// PBRT_CONSTEXPR int EndTransformBits = 1 << 1;
PBRT_CONSTEXPR int AllTransformsBits = (1 << MaxTransforms) - 1;
struct TransformSet {
    // TransformSet Public Methods
    Transform &operator[](int i) {
        CHECK_GE(i, 0);
        CHECK_LT(i, MaxTransforms);
        return t[i];
    }
    const Transform &operator[](int i) const {
        CHECK_GE(i, 0);
        CHECK_LT(i, MaxTransforms);
        return t[i];
    }
    friend TransformSet Inverse(const TransformSet &ts) {
        TransformSet tInv;
        for (int i = 0; i < MaxTransforms; ++i) tInv.t[i] = Inverse(ts.t[i]);
        return tInv;
    }
    bool IsAnimated() const {
        for (int i = 0; i < MaxTransforms - 1; ++i)
            if (t[i] != t[i + 1]) return true;
        return false;
    }

  private:
    Transform t[MaxTransforms];
};

struct RenderOptions {
    // RenderOptions Public Methods
    Integrator *MakeIntegrator() const;
    Scene *MakeScene();
    Camera *MakeCamera() const;

    // RenderOptions Public Data
    Float transformStartTime = 0, transformEndTime = 1;
    std::string FilterName = "box";
    ParamSet FilterParams;
    std::string FilmName = "image";
    ParamSet FilmParams;
    std::string SamplerName = "halton";
    ParamSet SamplerParams;
    std::string AcceleratorName = "bvh";
    ParamSet AcceleratorParams;
    std::string IntegratorName = "path";
    ParamSet IntegratorParams;
    std::string CameraName = "perspective";
    ParamSet CameraParams;
    TransformSet CameraToWorld;
    std::map<std::string, Medium*> namedMedia;
    std::vector<Light*> lights;
    std::vector<Primitive*> primitives;
    std::map<std::string, std::vector<Primitive*>> instances;
    std::vector<Primitive*> *currentInstance = nullptr;
    bool haveScatteringMedia = false;
};

// MaterialInstance represents both an instance of a material as well as
// the information required to create another instance of it (possibly with
// different parameters from the shape).
struct MaterialInstance {
    MaterialInstance() = default;
    MaterialInstance(const std::string &name, Material *mtl,
                     ParamSet params)
        : name(name), params(std::move(params)) 
    {
        material = mtl;
    }

    std::string name;
    Material* material;
    ParamSet params;
};

struct GraphicsState {
    // Graphics State Methods
    GraphicsState()
    // FIXME: would it be host?
    {
        floatTextures = new FloatTextureMap();
        spectrumTextures = new SpectrumTextureMap();
        namedMaterials = new NamedMaterialMap();
        ParamSet empty;
        TextureParams tp(empty, empty, *floatTextures, *spectrumTextures);
        Material* mtl(CreateMatteMaterial(tp));
        currentMaterial = new MaterialInstance("matte", mtl, ParamSet());
    }
    Material* GetMaterialForShape(const ParamSet &geomParams);
    MediumInterface CreateMediumInterface();

    // Graphics State
    std::string currentInsideMedium, currentOutsideMedium;

    // Updated after book publication: floatTextures, spectrumTextures, and
    // namedMaterials are all implemented using a "copy on write" approach
    // for more efficient GraphicsState management.  When state is pushed
    // in pbrtAttributeBegin(), we don't immediately make a copy of these
    // maps, but instead record that each one is shared.  Only if an item
    // is added to one is a unique copy actually made.
    using FloatTextureMap = std::map<std::string, std::shared_ptr<Texture<Float>>>;
    FloatTextureMap* floatTextures;
    bool floatTexturesShared = false;

    using SpectrumTextureMap = std::map<std::string, std::shared_ptr<Texture<Spectrum>>>;
    SpectrumTextureMap* spectrumTextures;
    bool spectrumTexturesShared = false;

    using NamedMaterialMap = std::map<std::string, std::shared_ptr<MaterialInstance>>;
    NamedMaterialMap* namedMaterials;
    bool namedMaterialsShared = false;

    MaterialInstance* currentMaterial;
    ParamSet areaLightParams;
    std::string areaLight;
    bool reverseOrientation = false;
};

STAT_MEMORY_COUNTER("Memory/TransformCache", transformCacheBytes);
STAT_PERCENT("Scene/TransformCache hits", nTransformCacheHits, nTransformCacheLookups);
STAT_INT_DISTRIBUTION("Scene/Probes per TransformCache lookup", transformCacheProbes);

// Note: TransformCache has been reimplemented and has a slightly different
// interface compared to the version described in the third edition of
// Physically Based Rendering.  The new version is more efficient in both
// space and memory, which is helpful for highly complex scenes.
//
// The new implementation uses a hash table to store Transforms (rather
// than a std::map, which generally uses a red-black tree).  Further,
// it doesn't always store the inverse of the transform; if a caller
// wants the inverse as well, they are responsible for storing it.
//
// The hash table size is always a power of two, allowing for the use of a
// bitwise AND to turn hash values into table offsets.  Quadratic probing
// is used when there is a hash collision.
class TransformCache {
  public:
    TransformCache()
        : hashTable(512), hashTableOccupancy(0) {}

    // TransformCache Public Methods
    Transform *Lookup(const Transform &t) {
        ++nTransformCacheLookups;

        int offset = Hash(t) & (hashTable.size() - 1);
        int step = 1;
        while (true) {
            // Keep looking until we find the Transform or determine that
            // it's not present.
            if (!hashTable[offset] || *hashTable[offset] == t)
                break;
            // Advance using quadratic probing.
            offset = (offset + step * step) & (hashTable.size() - 1);
            ++step;
        }
        ReportValue(transformCacheProbes, step);
        Transform *tCached = hashTable[offset];
        if (tCached)
            ++nTransformCacheHits;
        else {
            tCached = arena.Alloc<Transform>();
            *tCached = t;
            Insert(tCached);
        }
        return tCached;
    }

    void Clear() {
        transformCacheBytes += arena.TotalAllocated() + hashTable.size() * sizeof(Transform *);
        hashTable.clear();
        hashTable.resize(512);
        hashTableOccupancy = 0;
        arena.Reset();
    }

  private:
    void Insert(Transform *tNew);
    void Grow();

    static uint64_t Hash(const Transform &t) {
        const char *ptr = (const char *)(&t.GetMatrix());
        size_t size = sizeof(Matrix4x4);
        uint64_t hash = 14695981039346656037ull;
        while (size > 0) {
            hash ^= *ptr;
            hash *= 1099511628211ull;
            ++ptr;
            --size;
        }
        return hash;
    }

    // TransformCache Private Data
    std::vector<Transform *> hashTable;
    int hashTableOccupancy;
    MemoryArena arena;
};

void TransformCache::Insert(Transform *tNew) {
    if (++hashTableOccupancy == hashTable.size() / 2)
        Grow();

    int baseOffset = Hash(*tNew) & (hashTable.size() - 1);
    for (int nProbes = 0;; ++nProbes) {
        // Quadratic probing.
        int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
        if (hashTable[offset] == nullptr) {
            hashTable[offset] = tNew;
            return;
        }
    }
}

void TransformCache::Grow() {
    std::vector<Transform *> newTable(2 * hashTable.size());
    LOG(INFO) << "Growing transform cache hash table to " << newTable.size();

    // Insert current elements into newTable.
    for (Transform *tEntry : hashTable) {
        if (!tEntry) continue;

        int baseOffset = Hash(*tEntry) & (hashTable.size() - 1);
        for (int nProbes = 0;; ++nProbes) {
            // Quadratic probing.
            int offset = (baseOffset + nProbes/2 + nProbes*nProbes/2) & (hashTable.size() - 1);
            if (newTable[offset] == nullptr) {
                newTable[offset] = tEntry;
                break;
            }
        }
    }

    std::swap(hashTable, newTable);
}


// API Static Data
enum class APIState { Uninitialized, OptionsBlock, WorldBlock };
static APIState currentApiState = APIState::Uninitialized;
static TransformSet curTransform;
static uint32_t activeTransformBits = AllTransformsBits;
static std::map<std::string, TransformSet> namedCoordinateSystems;
static std::unique_ptr<RenderOptions> renderOptions;
static GraphicsState graphicsState;
static std::vector<GraphicsState> pushedGraphicsStates;
static std::vector<TransformSet> pushedTransforms;
static std::vector<uint32_t> pushedActiveTransformBits;
static TransformCache transformCache;
int catIndentCount = 0;

// API Forward Declarations
std::vector<Shape*> MakeShapes(const std::string &name,
                                               const Transform *ObjectToWorld,
                                               const Transform *WorldToObject,
                                               bool reverseOrientation,
                                               const ParamSet &paramSet);

// API Macros
#define VERIFY_INITIALIZED(func)                           \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::Uninitialized) {        \
        Error(                                             \
            "pbrtInit() must be before calling \"%s()\". " \
            "Ignoring.",                                   \
            func);                                         \
        return;                                            \
    } else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func)                             \
    VERIFY_INITIALIZED(func);                            \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&       \
        currentApiState == APIState::WorldBlock) {       \
        Error(                                           \
            "Options cannot be set inside world block; " \
            "\"%s\" not allowed.  Ignoring.",            \
            func);                                       \
        return;                                          \
    } else /* swallow trailing semicolon */
#define VERIFY_WORLD(func)                                   \
    VERIFY_INITIALIZED(func);                                \
    if (!(PbrtOptions.cat || PbrtOptions.toPly) &&           \
        currentApiState == APIState::OptionsBlock) {         \
        Error(                                               \
            "Scene description must be inside world block; " \
            "\"%s\" not allowed. Ignoring.",                 \
            func);                                           \
        return;                                              \
    } else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS(expr)           \
    for (int i = 0; i < MaxTransforms; ++i)   \
        if (activeTransformBits & (1 << i)) { \
            expr                              \
        }
#define WARN_IF_ANIMATED_TRANSFORM(func)                             \
    do {                                                             \
        if (curTransform.IsAnimated())                               \
            Warning(                                                 \
                "Animated transformations set; ignoring for \"%s\" " \
                "and using the start transform only",                \
                func);                                               \
    } while (false) /* swallow trailing semicolon */

// Object Creation Function Definitions
std::vector<Shape*> MakeShapes(const std::string &name,
                                               const Transform *object2world,
                                               const Transform *world2object,
                                               bool reverseOrientation,
                                               const ParamSet &paramSet) {
    std::vector<Shape*> shapes;
    Shape* s;
    if (name == "sphere")
        s = CreateSphereShape(object2world, world2object, reverseOrientation,
                              paramSet);
    if (s != nullptr) shapes.push_back(s);
    
    else
        Warning("Shape \"%s\" unknown.", name.c_str());
    return shapes;
}

STAT_COUNTER("Scene/Materials created", nMaterialsCreated);

Material* MakeMaterial(const std::string &name, const TextureParams &mp) {
    Material *material = nullptr;
    if (name == "" || name == "none")
        return nullptr;
    else if (name == "matte")
        material = CreateMatteMaterial(mp);
    else {
        Warning("Material \"%s\" unknown. Using \"matte\".", name.c_str());
        material = CreateMatteMaterial(mp);
    }

    if ((name == "subsurface" || name == "kdsubsurface") &&
        (renderOptions->IntegratorName != "path" &&
         (renderOptions->IntegratorName != "volpath")))
        Warning(
            "Subsurface scattering material \"%s\" used, but \"%s\" "
            "integrator doesn't support subsurface scattering. "
            "Use \"path\" or \"volpath\".",
            name.c_str(), renderOptions->IntegratorName.c_str());

    mp.ReportUnused();
    if (!material) Error("Unable to create material \"%s\"", name.c_str());
    else ++nMaterialsCreated;
    return material;
}

AreaLight* MakeAreaLight(const std::string &name,
                                         const Transform &light2world,
                                         const MediumInterface &mediumInterface,
                                         const ParamSet &paramSet,
                                         const Shape* shape) {
    AreaLight* area;
    if (name == "area" || name == "diffuse")
        area = CreateDiffuseAreaLight(light2world, mediumInterface.outside,
                                      paramSet, shape);
    else
        Warning("Area light \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return area;
}

Primitive* MakeAccelerator(
    const std::string &name,
    std::vector<Primitive*> prims,
    const ParamSet &paramSet) {
    Primitive* accel;
    if (name == "bvh")
        accel = CreateBVHAccelerator(std::move(prims), paramSet);
    else{
        Warning("Accelerator \"%s\" unknown.", name.c_str());
        exit(0);
    }
    paramSet.ReportUnused();
    return accel;
}

Camera *MakeCamera(const std::string &name, const ParamSet &paramSet,
                   const TransformSet &cam2worldSet, Float transformStart,
                   Float transformEnd, Film *film) {
    Camera *camera = nullptr;
    MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
    static_assert(MaxTransforms == 2,
                  "TransformCache assumes only two transforms");
    Transform *cam2world[2] = {
        transformCache.Lookup(cam2worldSet[0]),
        transformCache.Lookup(cam2worldSet[1])
    };
    AnimatedTransform animatedCam2World(cam2world[0], transformStart,
                                        cam2world[1], transformEnd);
    if (name == "perspective")
        camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film,
                                         mediumInterface.outside);
    else
        Warning("Camera \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return camera;
}

Sampler* MakeSampler(const std::string &name,
                                     const ParamSet &paramSet,
                                     const Film *film) {
    Sampler *sampler = nullptr;
    if (name == "halton")
        sampler = CreateHaltonSampler(paramSet, film->GetSampleBounds());
    else
        Warning("Sampler \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return sampler;
}

std::unique_ptr<Filter> MakeFilter(const std::string &name,
                                   const ParamSet &paramSet) {
    Filter *filter = nullptr;
    if (name == "box")
        filter = CreateBoxFilter(paramSet);
    else {
        Error("Filter \"%s\" unknown.", name.c_str());
        exit(1);
    }
    paramSet.ReportUnused();
    return std::unique_ptr<Filter>(filter);
}

Film *MakeFilm(const std::string &name, const ParamSet &paramSet,
               std::unique_ptr<Filter> filter) {
    Film *film = nullptr;
    if (name == "image")
        film = CreateFilm(paramSet, std::move(filter));
    else
        Warning("Film \"%s\" unknown.", name.c_str());
    paramSet.ReportUnused();
    return film;
}

// API Function Definitions
void pbrtInit(const Options &opt) {
    PbrtOptions = opt;
    // API Initialization
    if (currentApiState != APIState::Uninitialized)
        Error("pbrtInit() has already been called.");
    currentApiState = APIState::OptionsBlock;
    renderOptions.reset(new RenderOptions);
    graphicsState = GraphicsState();
    catIndentCount = 0;

    // General \pbrt Initialization
    ParallelInit();  // Threads must be launched before the profiler is
                     // initialized.
    InitProfiler();
}

void pbrtCleanup() {
    // API Cleanup
    if (currentApiState == APIState::Uninitialized)
        Error("pbrtCleanup() called without pbrtInit().");
    else if (currentApiState == APIState::WorldBlock)
        Error("pbrtCleanup() called while inside world block.");
    currentApiState = APIState::Uninitialized;
    ParallelCleanup();
    CleanupProfiler();
}

void pbrtIdentity() {
    VERIFY_INITIALIZED("Identity");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform();)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sIdentity\n", catIndentCount, "");
}

void pbrtTranslate(Float dx, Float dy, Float dz) {
    VERIFY_INITIALIZED("Translate");
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] *
                                            Translate(Vector3f(dx, dy, dz));)
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("%*sTranslate %.9g %.9g %.9g\n", catIndentCount, "", dx, dy,
               dz);
}

void pbrtLookAt(Float ex, Float ey, Float ez, Float lx, Float ly, Float lz,
                Float ux, Float uy, Float uz) {
    VERIFY_INITIALIZED("LookAt");
    Transform lookAt =
        LookAt(Point3f(ex, ey, ez), Point3f(lx, ly, lz), Vector3f(ux, uy, uz));
    FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookAt;);
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf(
            "%*sLookAt %.9g %.9g %.9g\n%*s%.9g %.9g %.9g\n"
            "%*s%.9g %.9g %.9g\n",
            catIndentCount, "", ex, ey, ez, catIndentCount + 8, "", lx, ly, lz,
            catIndentCount + 8, "", ux, uy, uz);
}

void pbrtFilm(const std::string &type, const ParamSet &params) {
    VERIFY_OPTIONS("Film");
    renderOptions->FilmParams = params;
    renderOptions->FilmName = type;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sFilm \"%s\" ", catIndentCount, "", type.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtSampler(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Sampler");
    renderOptions->SamplerName = name;
    renderOptions->SamplerParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sSampler \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtIntegrator(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Integrator");
    renderOptions->IntegratorName = name;
    renderOptions->IntegratorParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sIntegrator \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtCamera(const std::string &name, const ParamSet &params) {
    VERIFY_OPTIONS("Camera");
    renderOptions->CameraName = name;
    renderOptions->CameraParams = params;
    renderOptions->CameraToWorld = Inverse(curTransform);
    namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sCamera \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtWorldBegin() {
    VERIFY_OPTIONS("WorldBegin");
    currentApiState = APIState::WorldBlock;
    for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems["world"] = curTransform;
    if (PbrtOptions.cat || PbrtOptions.toPly)
        printf("\n\nWorldBegin\n\n");
}

void pbrtAttributeBegin() {
    VERIFY_WORLD("AttributeBegin");
    pushedGraphicsStates.push_back(graphicsState);
    graphicsState.floatTexturesShared = graphicsState.spectrumTexturesShared =
        graphicsState.namedMaterialsShared = true;
    pushedTransforms.push_back(curTransform);
    pushedActiveTransformBits.push_back(activeTransformBits);
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("\n%*sAttributeBegin\n", catIndentCount, "");
        catIndentCount += 4;
    }
}

void pbrtAttributeEnd() {
    VERIFY_WORLD("AttributeEnd");
    if (!pushedGraphicsStates.size()) {
        Error(
            "Unmatched pbrtAttributeEnd() encountered. "
            "Ignoring it.");
        return;
    }
    graphicsState = std::move(pushedGraphicsStates.back());
    pushedGraphicsStates.pop_back();
    curTransform = pushedTransforms.back();
    pushedTransforms.pop_back();
    activeTransformBits = pushedActiveTransformBits.back();
    pushedActiveTransformBits.pop_back();
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        catIndentCount -= 4;
        printf("%*sAttributeEnd\n", catIndentCount, "");
    }
}

void pbrtMaterial(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Material");
    ParamSet emptyParams;
    TextureParams mp(params, emptyParams, *graphicsState.floatTextures,
                     *graphicsState.spectrumTextures);
    Material* mtl = MakeMaterial(name, mp);
    graphicsState.currentMaterial = new MaterialInstance(name, mtl, params);

    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sMaterial \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtAreaLightSource(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("AreaLightSource");
    graphicsState.areaLight = name;
    graphicsState.areaLightParams = params;
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sAreaLightSource \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }
}

void pbrtShape(const std::string &name, const ParamSet &params) {
    VERIFY_WORLD("Shape");
    std::vector<Primitive*> prims;
    std::vector<AreaLight*> areaLights;
    if (PbrtOptions.cat || (PbrtOptions.toPly && name != "trianglemesh")) {
        printf("%*sShape \"%s\" ", catIndentCount, "", name.c_str());
        params.Print(catIndentCount);
        printf("\n");
    }

    if (!curTransform.IsAnimated()) {
        // Initialize _prims_ and _areaLights_ for static shape

        // Create shapes for shape _name_
        Transform *ObjToWorld = transformCache.Lookup(curTransform[0]);
        Transform *WorldToObj = transformCache.Lookup(Inverse(curTransform[0]));
        std::vector<Shape*> shapes =
            MakeShapes(name, ObjToWorld, WorldToObj,
                       graphicsState.reverseOrientation, params);
        if (shapes.empty()) return;
        Material* mtl = graphicsState.GetMaterialForShape(params);
        params.ReportUnused();
        MediumInterface mi = graphicsState.CreateMediumInterface();
        prims.reserve(shapes.size());
        GeometricPrimitive* ptr;
        cudaMallocManaged(&ptr, sizeof(GeometricPrimitive) * shapes.size());
        for (auto s : shapes) {
            // Possibly create area light for shape
            AreaLight* area;
            if (graphicsState.areaLight != "") {

                area = MakeAreaLight(graphicsState.areaLight, curTransform[0],
                                     mi, graphicsState.areaLightParams, s);
                if (area) areaLights.push_back(area);
            }
            prims.push_back(
                new(ptr) GeometricPrimitive(s, mtl, area, mi));
            ptr++;
        }
    } else {
        printf("Animated Scene not supported.\n");
        assert(false);
    }
    // Add _prims_ and _areaLights_ to scene or current instance
    if (renderOptions->currentInstance) {
        if (areaLights.size())
            Warning("Area lights not supported with object instancing");
        renderOptions->currentInstance->insert(
            renderOptions->currentInstance->end(), prims.begin(), prims.end());
    } else {
        renderOptions->primitives.insert(renderOptions->primitives.end(),
                                         prims.begin(), prims.end());
        if (areaLights.size())
            renderOptions->lights.insert(renderOptions->lights.end(),
                                         areaLights.begin(), areaLights.end());
    }
}

// Attempt to determine if the ParamSet for a shape may provide a value for
// its material's parameters. Unfortunately, materials don't provide an
// explicit representation of their parameters that we can query and
// cross-reference with the parameter values available from the shape.
//
// Therefore, we'll apply some "heuristics".
bool shapeMaySetMaterialParameters(const ParamSet &ps) {
    for (const auto &param : ps.textures)
        // Any texture other than one for an alpha mask is almost certainly
        // for a Material (or is unused!).
        if (param->name != "alpha" && param->name != "shadowalpha")
            return true;

    // Special case spheres, which are the most common non-mesh primitive.
    for (const auto &param : ps.floats)
        if (param->nValues == 1 && param->name != "radius")
            return true;

    // Extra special case strings, since plymesh uses "filename", curve "type",
    // and loopsubdiv "scheme".
    for (const auto &param : ps.strings)
        if (param->nValues == 1 && param->name != "filename" &&
            param->name != "type" && param->name != "scheme")
            return true;

    // For all other parameter types, if there is a single value of the
    // parameter, assume it may be for the material. This should be valid
    // (if conservative), since no materials currently take array
    // parameters.
    for (const auto &param : ps.bools)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.ints)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector2fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.point3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.vector3fs)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.normals)
        if (param->nValues == 1)
            return true;
    for (const auto &param : ps.spectra)
        if (param->nValues == 1)
            return true;

    return false;
}

Material* GraphicsState::GetMaterialForShape(
    const ParamSet &shapeParams) {
    CHECK(currentMaterial);
    if (shapeMaySetMaterialParameters(shapeParams)) {
        // Only create a unique material for the shape if the shape's
        // parameters are (apparently) going to provide values for some of
        // the material parameters.
        TextureParams mp(shapeParams, currentMaterial->params, *floatTextures,
                         *spectrumTextures);
        return MakeMaterial(currentMaterial->name, mp);
    } else
        return currentMaterial->material;
}

MediumInterface GraphicsState::CreateMediumInterface() {
    MediumInterface m;
    if (currentInsideMedium != "") {
        if (renderOptions->namedMedia.find(currentInsideMedium) !=
            renderOptions->namedMedia.end())
            m.inside = renderOptions->namedMedia[currentInsideMedium];
        else
            Error("Named medium \"%s\" undefined.",
                  currentInsideMedium.c_str());
    }
    if (currentOutsideMedium != "") {
        if (renderOptions->namedMedia.find(currentOutsideMedium) !=
            renderOptions->namedMedia.end())
            m.outside = renderOptions->namedMedia[currentOutsideMedium];
        else
            Error("Named medium \"%s\" undefined.",
                  currentOutsideMedium.c_str());
    }
    return m;
}

STAT_COUNTER("Scene/Object instances created", nObjectInstancesCreated);

STAT_COUNTER("Scene/Object instances used", nObjectInstancesUsed);

void pbrtWorldEnd() {
    VERIFY_WORLD("WorldEnd");
    // Ensure there are no pushed graphics states
    while (pushedGraphicsStates.size()) {
        Warning("Missing end to pbrtAttributeBegin()");
        pushedGraphicsStates.pop_back();
        pushedTransforms.pop_back();
    }
    while (pushedTransforms.size()) {
        Warning("Missing end to pbrtTransformBegin()");
        pushedTransforms.pop_back();
    }

    // Create scene and render
    if (PbrtOptions.cat || PbrtOptions.toPly) {
        printf("%*sWorldEnd\n", catIndentCount, "");
    } else {
        std::unique_ptr<Integrator> integrator(renderOptions->MakeIntegrator());
        std::unique_ptr<Scene> scene(renderOptions->MakeScene());

        // This is kind of ugly; we directly override the current profiler
        // state to switch from parsing/scene construction related stuff to
        // rendering stuff and then switch it back below. The underlying
        // issue is that all the rest of the profiling system assumes
        // hierarchical inheritance of profiling state; this is the only
        // place where that isn't the case.
        CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::SceneConstruction));
        ProfilerState = ProfToBits(Prof::IntegratorRender);

        if (scene && integrator) integrator->Render(*scene);

        CHECK_EQ(CurrentProfilerState(), ProfToBits(Prof::IntegratorRender));
        ProfilerState = ProfToBits(Prof::SceneConstruction);
    }

    // Clean up after rendering. Do this before reporting stats so that
    // destructors can run and update stats as needed.
    graphicsState = GraphicsState();
    transformCache.Clear();
    currentApiState = APIState::OptionsBlock;
    //ImageTexture<Float, Float>::ClearCache();
    //ImageTexture<RGBSpectrum, Spectrum>::ClearCache();
    renderOptions.reset(new RenderOptions);

    if (!PbrtOptions.cat && !PbrtOptions.toPly) {
        MergeWorkerThreadStats();
        ReportThreadStats();
        if (!PbrtOptions.quiet) {
            PrintStats(stdout);
            ReportProfilerResults(stdout);
            ClearStats();
            ClearProfiler();
        }
    }

    for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    activeTransformBits = AllTransformsBits;
    namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
                                 namedCoordinateSystems.end());
}

Scene *RenderOptions::MakeScene() {
    Primitive* accelerator =
        MakeAccelerator(AcceleratorName, std::move(primitives), AcceleratorParams);
    // if (!accelerator) accelerator = std::make_shared<BVHAccel>(primitives);
    Scene *scene = new Scene(accelerator, lights.data());
    // Erase primitives and lights from _RenderOptions_
    primitives.clear();
    lights.clear();
    return scene;
}

Integrator *RenderOptions::MakeIntegrator() const {
    const Camera* camera(MakeCamera());
    if (!camera) {
        Error("Unable to create camera");
        return nullptr;
    }

    Sampler* sampler =
        MakeSampler(SamplerName, SamplerParams, camera->film);
    if (!sampler) {
        Error("Unable to create sampler.");
        return nullptr;
    }

    Integrator *integrator = nullptr;
    if (IntegratorName == "path")
        integrator = CreatePathIntegrator(IntegratorParams, sampler, camera);
    else {
        Error("Integrator \"%s\" unknown.", IntegratorName.c_str());
        return nullptr;
    }

    if (renderOptions->haveScatteringMedia && IntegratorName != "volpath" &&
        IntegratorName != "bdpt" && IntegratorName != "mlt") {
        Warning(
            "Scene has scattering media but \"%s\" integrator doesn't support "
            "volume scattering. Consider using \"volpath\", \"bdpt\", or "
            "\"mlt\".", IntegratorName.c_str());
    }

    IntegratorParams.ReportUnused();
    // Warn if no light sources are defined
    if (lights.empty())
        Warning(
            "No light sources defined in scene; "
            "rendering a black image.");
    return integrator;
}

Camera *RenderOptions::MakeCamera() const {
    std::unique_ptr<Filter> filter = MakeFilter(FilterName, FilterParams);
    Film *film = MakeFilm(FilmName, FilmParams, std::move(filter));
    if (!film) {
        Error("Unable to create film.");
        return nullptr;
    }
    Camera *camera = pbrt::MakeCamera(CameraName, CameraParams, CameraToWorld,
                                  renderOptions->transformStartTime,
                                  renderOptions->transformEndTime, film);
    return camera;
}

}  // namespace pbrt
