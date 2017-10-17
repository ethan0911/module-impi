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

/*! ISCP side code of the impi module: this code doesn't know anything
    specific about the underlying "thing" that defines voxels/voxels -
    it only gets a list of voxel "references" (each on ebeing a 64-bit
    integer), and call-back function (to the C side) to query the
    bounds and/or voxel values of said voxels. 

    All the actual implemntation of what a voxel looks like - what
    format or voxel type, where they come from (explicitly from a list
    of voxels extracted with another tool, or implicitly in a 3D
    strcutured volume, etc), lies on the C side - this code ONLY
    builds a bvh over the voxel refs, and intersects them as required
    during traversal */


/*! function to get the 3D world bounds of a given voxel reference, to be computed on the C++ side */
unmasked extern "C" externC_getVoxelBounds(uniform box3fa &voxelBounds,
                                           void *uniform c_self,
                                           const uniform uint64 voxelRef);

/*! function to get the corner voxel values of a given voxel (referend
    with given voxelRef), to be computed on the C++ side */
unmasked extern "C" externC_getVoxel(uniform Voxel &voxel,
                                     void *uniform c_self,
                                     const uniform uint64 voxelRef);


struct Impi {
  /*! inherit from "Geometry" class: since ISPC doesn't support
      inheritance we simply put the "parent" class as the first
      member; this way any typecast to the parent class will get the
      right members (including 'virtual' function pointers, etc) */
  Geometry   super;

  float isoValue;
  float *voxelArray;
  
  /*! for the case where we build an embree bvh over the hot voxels,
      this is the list of all voxels that are hot (each one is one prim
      in the embree BVH */
  uint64 *uniform activeVoxelRefs;

  /*! for the embree bvh over active voxels case this is the c-handle to
      the c-side volume. this is _probably_ a C++-side virtual class
      that implements getvoxelbounds and getvoxel */
  void *uniform c_self;

  /*! todo - add getVoxel and getVoxelBounds as member function pointers
      (and let c++ side pass them on constructor), rather than as
      global functions */
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
unmasked void Impi_bounds_testVoxel(Impi *uniform self,
                                   uniform size_t primID,
                                   uniform box3fa &bbox)
{
  bbox.lower = make_vec3fa(0.f);
  bbox.upper = make_vec3fa(1.f);
}



/*! embree callback to compute voxel bounds - we simply pass that on to
    the C side (yes, might eventually do this completely on the
    C-side, without first going through this callback, too.... but
    don't (yet) want to deal with having to have some callbacks like
    on the C side (getbounds) and others (intersectVoxel) on the ISPC
    side) */
unmasked void Impi_bounds(Impi *uniform self,
                          uniform size_t primID,
                          uniform box3fa &bbox)
{
  externC_getVoxelBounds(bbox,self->c_self,self->activeVoxelRefs[primID]);
}



/*! intersect callback for the 'embree over hot voxels' method */
void Impi_intersect(Impi *uniform self,
                             varying Ray &ray,
                             uniform size_t primID)
{
  uniform Voxel  voxel;
  externC_getVoxel(voxel,self->c_self,self->activeVoxelRefs[primID]);
  if (actualVoxelIntersect(ray,voxel,self->isoValue)) {
    ray.geomID = self->super.geomID;
    ray.primID = primID;
  }
  return;
}


export void Impi_finalize(void   *uniform _self,
                          void   *uniform _model,
                          uint64 *uniform activeVoxelRefs,
                          uint64  uniform numActiveVoxelRefs,
                          void   *uniform c_self,
                          uniform float   isoValue)
{
  // first, typecast to our 'real' type. since ispc can't export real
  // types to c we have to pass 'self' in as a void*, and typecast
  // here.
  Impi *uniform self = (Impi *uniform)_self;
  
  // set our internal data.
  self->isoValue   = isoValue;
  self->activeVoxelRefs = activeVoxelRefs;
  self->c_self      = c_self;
  
  // ... and let embree build a bvh, with 'numPatches' primitmives and
  // our function pointers for ray-prim interseciton and
  // prim-boundingbox computation callsbacks
  Model *uniform model  = (Model *uniform)_model;
  
  // create a new embree geometry with numpathces prims, in the model
  // that this goemetry is in.
  uint32 uniform geomID = rtcNewUserGeometry(model->embreeSceneHandle,numActiveVoxelRefs);
  
  // set 'us' as user data (this will be the first arg in intersect()
  // and computebounds() callbacks
  rtcSetUserData(model->embreeSceneHandle,geomID,self);
  
  // set function pointers:
  rtcSetBoundsFunction(model->embreeSceneHandle,geomID,
                       (uniform RTCBoundsFunc)&Impi_bounds);
  rtcSetIntersectFunction(model->embreeSceneHandle,geomID,
                          (uniform RTCIntersectFuncVarying)&Impi_intersect);
  rtcSetOccludedFunction(model->embreeSceneHandle,geomID,
                         (uniform RTCOccludedFuncVarying)&Impi_intersect);
}
                                    

