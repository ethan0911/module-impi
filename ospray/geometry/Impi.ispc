// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

// ospray
#include "math/vec.ih"
#include "math/box.ih"
#include "math/LinearSpace.ih"
#include "common/Ray.ih"
#include "common/Model.ih"
#include "ospray/geometry/Geometry.ih"
// embree
#include "embree2/rtcore.isph"
#include "embree2/rtcore_scene.isph"
#include "embree2/rtcore_geometry_user.isph"

#include "Voxel.ih"

struct Impi {
  /*! inherit from "Geometry" class: since ISPC doesn't support
      inheritance we simply put the "parent" class as the first
      member; this way any typecast to the parent class will get the
      right members (including 'virtual' function pointers, etc) */
  Geometry   super;

  float isoValue;
  float *voxelArray;
  vec3f rcpDims; // rcp(dims)

  /*! for the case where we build an embree bvh over the hot cells,
      this is the list of all cells that are hot (each one is one prim
      in the embree BVH */
  uint64 *uniform hotCellRefs;

  /*! for the embree bvh over hot cells case this is the c-handle to
      the c-side volume */
  void *uniform volume;
};

static void Impi_postIntersect(uniform Geometry *uniform geometry,
                               uniform Model *uniform model,
                               varying DifferentialGeometry &dg,
                               const varying Ray &ray,
                               uniform int64 flags)
{
  Impi *uniform self = (Impi *uniform)geometry;
  dg.Ng = dg.Ns = normalize(ray.Ng);

  if (flags & DG_MATERIALID) {
    dg.materialID = 0;
    dg.material   = self->super.material;
  }
}


/*! create a new ispc equivalent for the C++ class (cppEquivalent
    points back, if we ever need this. this returns a void * that the
    c++ side can use to call back to us when needed. Note this return
    values _has_ to be a untyped 'void*' because ISPC cannot currently
    export 'real' class types - with function pointers etcpp - to a c
    header file */
export void *uniform Impi_create(void *uniform cppEquivalent)
{
  Impi *uniform self = uniform new uniform Impi;
  
  /* recursively call the 'inherited' constructor, and pass our
     'overridden' functions for intersect and postintersect. note we
     do not have to set the 'intersect' and 'occluded' functions
     because we use embree to do the traversal for us (see
     Impi_set()) below */
  Geometry_Constructor(&self->super,cppEquivalent,
                       Impi_postIntersect,
                       NULL,0,NULL);
  return self;
}

export void Impi_destroy(void *uniform _self)
{
  /* _actually_ this should also destroy the created embree geometry
     here; not doing this to keep the code as small as possible, but
     usually this is where it should be done */
  Impi *uniform self = (Impi *uniform)_self;
  delete self;
}


/*! to be able to have embree build a data structure over us we have
    to be able to tell it the bounding box for each primitimve. this
    is the callback embree calls for querying the (uniform) bounding
    box for a single given primitmive 'primID' */
unmasked void Impi_bounds_testCell(Impi *uniform self,
                                   uniform size_t primID,
                                   uniform box3fa &bbox)
{
  bbox.lower = make_vec3fa(0.f);
  bbox.upper = make_vec3fa(1.f);
}




#define CELL_REF_BITS 21
#define CELL_REF_MASK ((1<<CELL_REF_BITS)-1)

inline uniform vec3i extractCellRef(const uniform uint64 cellRef)
{
  return make_vec3i((int)(cellRef & CELL_REF_MASK),
                    (int)((cellRef >> CELL_REF_BITS) & CELL_REF_MASK),
                    (int)((cellRef >> (2*CELL_REF_BITS)) & CELL_REF_MASK));
}

unmasked void Impi_bounds_hotCells(Impi *uniform self,
                                   uniform size_t primID,
                                   uniform box3fa &bbox)
{
  uniform vec3i cellIdx = extractCellRef(self->hotCellRefs[primID]);
  bbox.lower = make_vec3fa(make_vec3f(cellIdx)*self->rcpDims);
  bbox.upper = bbox.lower + make_vec3fa(self->rcpDims);
}





/*! intersect callback for the simple test-cell example */
void Impi_intersect_testCell(Impi *uniform self,
                             varying Ray &ray,
                             uniform size_t primID)
{
  uniform Voxel *uniform v = (uniform Voxel *uniform)self->voxelArray;
  if (unitVoxelIntersect(ray,*v,self->isoValue)) {
    ray.geomID = self->super.geomID;
    ray.primID = primID;
    const vec3f P = getPoint(ray,ray.t);
    // print("hit\n %\n %\n %\n",P.x,P.y,P.z);
    ray.Ng = normalize(gradient(*v,P));
    // print("grad\n %\n %\n %\n",ray.Ng.x,ray.Ng.y,ray.Ng.z);
    // ray.Ng = make_vec3f(1,0,0);
  }
}


unmasked extern "C" externC_getCell(uniform Voxel &cell,
                                     void *uniform volume,
                                     const uniform vec3i &cellIdx);


/*! intersect callback for the 'embree over hot cells' method */
void Impi_intersect_hotCells(Impi *uniform self,
                             varying Ray &ray,
                             uniform size_t primID)
{
  uniform Voxel voxel;
  const uniform vec3i cellIdx = extractCellRef(self->hotCellRefs[primID]);
  externC_getCell(voxel,self->volume,cellIdx);
  if (actualVoxelIntersect(ray,voxel,cellIdx,self->rcpDims,self->isoValue)) {
    ray.geomID = self->super.geomID;
    ray.primID = primID;
  }
  return;
  // uniform Voxel *uniform v = (uniform Voxel *uniform)self->voxelArray;
  // if (unitVoxelIntersect(ray,*v,self->isoValue)) {
  //   ray.geomID = self->super.geomID;
  //   ray.primID = primID;
  //   const vec3f P = getPoint(ray,ray.t);
  //   // print("hit\n %\n %\n %\n",P.x,P.y,P.z);
  //   ray.Ng = normalize(gradient(*v,P));
  //   // print("grad\n %\n %\n %\n",ray.Ng.x,ray.Ng.y,ray.Ng.z);
  //   // ray.Ng = make_vec3f(1,0,0);
  // }
}



/*! the function called by Impi::commit() to pass the
  actual user data over. for simplicity we assume that this functoin
  gets called exactly once; in practice a user may also _change_ an
  existing class'es data members, but handling this would make this
  code significantly larger */
export void Impi_finalize_testCell(void *uniform _self,
                                   void *uniform _model,
                                   float *uniform voxelArray,
                                   uniform float isoValue)
{
  // first, typecast to our 'real' type. since ispc can't export real
  // types to c we have to pass 'self' in as a void*, and typecast
  // here.
  Impi *uniform self = (Impi *uniform)_self;
  
  // set our internal data.
  self->voxelArray = voxelArray;
  self->isoValue   = isoValue;
  
  // ... and let embree build a bvh, with 'numPatches' primitmives and
  // our function pointers for ray-prim interseciton and
  // prim-boundingbox computation callsbacks
  Model *uniform model  = (Model *uniform)_model;

  // create a new embree geometry with numpathces prims, in the model
  // that this goemetry is in.
  uint32 uniform geomID = rtcNewUserGeometry(model->embreeSceneHandle,1);
  
  // set 'us' as user data (this will be the first arg in intersect()
  // and computebounds() callbacks
  rtcSetUserData(model->embreeSceneHandle,geomID,self);
  
  // set function pointers:
  rtcSetBoundsFunction(model->embreeSceneHandle,geomID,
                       (uniform RTCBoundsFunc)&Impi_bounds_testCell);
  rtcSetIntersectFunction(model->embreeSceneHandle,geomID,
                          (uniform RTCIntersectFuncVarying)&Impi_intersect_testCell);
  rtcSetOccludedFunction(model->embreeSceneHandle,geomID,
                         (uniform RTCOccludedFuncVarying)&Impi_intersect_testCell);
}
                                    

export void Impi_finalize_embreeBVHoverHotCells(void *uniform _self,
                                                void *uniform _model,
                                                uint64 *uniform hotCellRefs,
                                                uint64 uniform numHotCellRefs,
                                                uniform vec3i &dims,
                                                void *uniform _volume,
                                                uniform float isoValue)
{
  // first, typecast to our 'real' type. since ispc can't export real
  // types to c we have to pass 'self' in as a void*, and typecast
  // here.
  Impi *uniform self = (Impi *uniform)_self;
  
  // set our internal data.
  self->isoValue   = isoValue;
  self->rcpDims  = rcp(make_vec3f(dims)-make_vec3f(1.f));
  self->hotCellRefs = hotCellRefs;
  self->volume      = _volume;
  
  // ... and let embree build a bvh, with 'numPatches' primitmives and
  // our function pointers for ray-prim interseciton and
  // prim-boundingbox computation callsbacks
  Model *uniform model  = (Model *uniform)_model;
  
  // create a new embree geometry with numpathces prims, in the model
  // that this goemetry is in.
  uint32 uniform geomID = rtcNewUserGeometry(model->embreeSceneHandle,numHotCellRefs);
  
  // set 'us' as user data (this will be the first arg in intersect()
  // and computebounds() callbacks
  rtcSetUserData(model->embreeSceneHandle,geomID,self);
  
  // set function pointers:
  rtcSetBoundsFunction(model->embreeSceneHandle,geomID,
                       (uniform RTCBoundsFunc)&Impi_bounds_hotCells);
  rtcSetIntersectFunction(model->embreeSceneHandle,geomID,
                          (uniform RTCIntersectFuncVarying)&Impi_intersect_hotCells);
  rtcSetOccludedFunction(model->embreeSceneHandle,geomID,
                         (uniform RTCOccludedFuncVarying)&Impi_intersect_hotCells);
}
                                    

