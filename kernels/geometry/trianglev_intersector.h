// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "triangle.h"
#include "triangle_intersector_pluecker.h"
#include "triangle_intersector_moeller.h"
#include "triangle_intersector_woop.h"

namespace embree
{
  namespace isa
  {
    template<int M>
    struct MoellerRefinedHitM
    {
      __forceinline MoellerRefinedHitM(const vbool<M>& valid, const Vec3vf<M>& v0, const Vec3vf<M>& v1, const Vec3vf<M>& v2, const Vec3vf<M>& D)
        : valid(valid), v0(v0), v1(v1), v2(v2), D(D) {}

      __forceinline void finalize()
      {
        const Vec3vf<M> e0 = v2-v0;
        const Vec3vf<M> e1 = v0-v1;
        const Vec3vf<M> e2 = v1-v2;
        const vfloat<M> U = dot(cross(e0,v2+v0),D);
        const vfloat<M> V = dot(cross(e1,v0+v1),D);
        const vfloat<M> W = dot(cross(e2,v1+v2),D);
        const vfloat<M> UVW = U+V+W;
        vNg = stable_triangle_normal(e0,e1,e2);
        const vfloat<M> den = twice(dot(vNg,D));
        const vfloat<M> T = twice(dot(v0,vNg));
        vt = rcp(den)*T;
        const vbool<M> invalid = abs(UVW) < min_rcp_input;
        const vfloat<M> rcpUVW = select(invalid,vfloat<M>(0.0f),rcp(UVW));
        vu = min(U * rcpUVW,1.0f);
        vv = min(V * rcpUVW,1.0f);
      }

      __forceinline Vec2vf<M> uv() const { return Vec2vf<M>(vu,vv); }
      __forceinline vfloat<M> t () const { return vt; }
      __forceinline Vec3vf<M> Ng() const { return vNg; }

      __forceinline Vec2f uv (const size_t i) const { return Vec2f(vu[i],vv[i]); }
      __forceinline float t  (const size_t i) const { return vt[i]; }
      __forceinline Vec3fa Ng(const size_t i) const { return Vec3fa(vNg.x[i],vNg.y[i],vNg.z[i]); }

    public:
      vbool<M> valid;
      Vec3vf<M> v0, v1, v2, D;
      vfloat<M> vu, vv, vt;
      Vec3vf<M> vNg;
    };

    template<int M>
    struct MoellerRefinedIntersector1
    {
      __forceinline MoellerRefinedIntersector1() {}
      __forceinline MoellerRefinedIntersector1(const Ray& ray, const void* ptr) {}

      template<typename Epilog>
      __forceinline bool intersect(Ray& ray, const Vec3vf<M>& tri_v0, const Vec3vf<M>& tri_v1, const Vec3vf<M>& tri_v2, const Epilog& epilog) const
      {
        auto mapUV = UVIdentity<M>();
        MoellerTrumboreHitM<M,UVIdentity<M>> mhit(mapUV);
        MoellerTrumboreIntersector1<M> moeller;
        if (likely(moeller.intersect(ray,tri_v0,tri_v1,tri_v2,mapUV,mhit)))
        {
          const Vec3vf<M> O = Vec3vf<M>((Vec3fa)ray.org);
          const Vec3vf<M> D = Vec3vf<M>((Vec3fa)ray.dir);
          MoellerRefinedHitM<M> hit(mhit.valid,tri_v0-O,tri_v1-O,tri_v2-O,D);
          return epilog(hit.valid,hit);
        }
        return false;
      }

      template<typename Epilog>
      __forceinline bool occluded(Ray& ray, const Vec3vf<M>& tri_v0, const Vec3vf<M>& tri_v1, const Vec3vf<M>& tri_v2, const Epilog& epilog) const
      {
        MoellerTrumboreIntersector1<M> moeller;
        return moeller.intersect(ray,tri_v0,tri_v1,tri_v2,epilog);
      }
    };

    template<int M, bool filter>
    struct TriangleMvIntersector1MoellerRefined
    {
      typedef TriangleMv<M> Primitive;
      typedef MoellerRefinedIntersector1<M> Precalculations;

      static __forceinline void intersect(Precalculations& pre, RayHit& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        pre.intersect(ray,tri.v0,tri.v1,tri.v2,Intersect1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }

      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return pre.occluded(ray,tri.v0,tri.v1,tri.v2,Occluded1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }

      static __forceinline bool pointQuery(PointQuery* query, PointQueryContext* context, const Primitive& tri)
      {
        return PrimitivePointQuery1<Primitive>::pointQuery(query, context, tri);
      }
    };


    /*! Intersects M triangles with 1 ray */
    template<int M, bool filter>
    struct TriangleMvIntersector1Moeller
    {
      typedef TriangleMv<M> Primitive;
      typedef MoellerTrumboreIntersector1<M> Precalculations;

      /*! Intersect a ray with M triangles and updates the hit. */
      static __forceinline void intersect(Precalculations& pre, RayHit& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        pre.intersect(ray,tri.v0,tri.v1,tri.v2,/*UVIdentity<M>(),*/Intersect1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }

      /*! Test if the ray is occluded by one of the M triangles. */
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return pre.intersect(ray,tri.v0,tri.v1,tri.v2,/*UVIdentity<M>(),*/Occluded1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }
      
      static __forceinline bool pointQuery(PointQuery* query, PointQueryContext* context, const Primitive& tri)
      {
        return PrimitivePointQuery1<Primitive>::pointQuery(query, context, tri);
      }
    };


    template<int M, bool filter>
    struct TriangleMvIntersector1Woop
    {
      typedef TriangleMv<M> Primitive;
      typedef WoopIntersector1<M> intersec;
      typedef WoopPrecalculations1<M> Precalculations;

      /*! Intersect a ray with M triangles and updates the hit. */
      static __forceinline void intersect(const Precalculations& pre, RayHit& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        intersec::intersect(ray,pre,tri.v0,tri.v1,tri.v2,Intersect1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }

      /*! Test if the ray is occluded by one of the M triangles. */
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return intersec::intersect(ray,pre,tri.v0,tri.v1,tri.v2,Occluded1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }
      
      static __forceinline bool pointQuery(PointQuery* query, PointQueryContext* context, const Primitive& tri)
      {
        return PrimitivePointQuery1<Primitive>::pointQuery(query, context, tri);
      }
    };


    /*! Intersects M triangles with K rays */
    template<int M, int K, bool filter>
    struct TriangleMvIntersectorKMoeller
    {
      typedef TriangleMv<M> Primitive;
      typedef MoellerTrumboreIntersectorK<M,K> Precalculations;

      /*! Intersects K rays with M triangles. */
      static __forceinline void intersect(const vbool<K>& valid_i, Precalculations& pre, RayHitK<K>& ray, RayQueryContext* context, const Primitive& tri)
      {
        for (size_t i=0; i<M; i++)
        {
          if (!tri.valid(i)) break;
          STAT3(normal.trav_prims,1,popcnt(valid_i),K);
          const Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0,i);
          const Vec3vf<K> v1 = broadcast<vfloat<K>>(tri.v1,i);
          const Vec3vf<K> v2 = broadcast<vfloat<K>>(tri.v2,i);
          pre.intersectK(valid_i,ray,v0,v1,v2,/*UVIdentity<K>(),*/IntersectKEpilogM<M,K,filter>(ray,context,tri.geomID(),tri.primID(),i));
        }
      }

      /*! Test for K rays if they are occluded by any of the M triangles. */
      static __forceinline vbool<K> occluded(const vbool<K>& valid_i, Precalculations& pre, RayK<K>& ray, RayQueryContext* context, const Primitive& tri)
      {
        vbool<K> valid0 = valid_i;

        for (size_t i=0; i<M; i++)
        {
          if (!tri.valid(i)) break;
          STAT3(shadow.trav_prims,1,popcnt(valid_i),K);
          const Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0,i);
          const Vec3vf<K> v1 = broadcast<vfloat<K>>(tri.v1,i);
          const Vec3vf<K> v2 = broadcast<vfloat<K>>(tri.v2,i);
          pre.intersectK(valid0,ray,v0,v1,v2,/*UVIdentity<K>(),*/OccludedKEpilogM<M,K,filter>(valid0,ray,context,tri.geomID(),tri.primID(),i));
          if (none(valid0)) break;
        }
        return !valid0;
      }
      
      /*! Intersect a ray with M triangles and updates the hit. */
      static __forceinline void intersect(Precalculations& pre, RayHitK<K>& ray, size_t k, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        pre.intersect(ray,k,tri.v0,tri.v1,tri.v2,/*UVIdentity<M>(),*/Intersect1KEpilogM<M,K,filter>(ray,k,context,tri.geomID(),tri.primID())); //FIXME: M
      }

      /*! Test if the ray is occluded by one of the M triangles. */
      static __forceinline bool occluded(Precalculations& pre, RayK<K>& ray, size_t k, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return pre.intersect(ray,k,tri.v0,tri.v1,tri.v2,/*UVIdentity<M>(),*/Occluded1KEpilogM<M,K,filter>(ray,k,context,tri.geomID(),tri.primID())); //FIXME: M
      }
    };

    /*! Intersects M triangles with 1 ray */
    template<int M, bool filter>
    struct TriangleMvIntersector1Pluecker
    {
      typedef TriangleMv<M> Primitive;
      typedef PlueckerIntersector1<M> Precalculations;

      /*! Intersect a ray with M triangles and updates the hit. */
      static __forceinline void intersect(Precalculations& pre, RayHit& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        pre.intersect(ray,tri.v0,tri.v1,tri.v2,UVIdentity<M>(),Intersect1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }

      /*! Test if the ray is occluded by one of the M triangles. */
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return pre.intersect(ray,tri.v0,tri.v1,tri.v2,UVIdentity<M>(),Occluded1EpilogM<M,filter>(ray,context,tri.geomID(),tri.primID()));
      }
      
      static __forceinline bool pointQuery(PointQuery* query, PointQueryContext* context, const Primitive& tri)
      {
        return PrimitivePointQuery1<Primitive>::pointQuery(query, context, tri);
      }
    };

    /*! Intersects M triangles with K rays */
    template<int M, int K, bool filter>
    struct TriangleMvIntersectorKPluecker
    {
      typedef TriangleMv<M> Primitive;
      typedef PlueckerIntersectorK<M,K> Precalculations;

      /*! Intersects K rays with M triangles. */
      static __forceinline void intersect(const vbool<K>& valid_i, Precalculations& pre, RayHitK<K>& ray, RayQueryContext* context, const Primitive& tri)
      {
        for (size_t i=0; i<M; i++)
        {
          if (!tri.valid(i)) break;
          STAT3(normal.trav_prims,1,popcnt(valid_i),K);
          const Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0,i);
          const Vec3vf<K> v1 = broadcast<vfloat<K>>(tri.v1,i);
          const Vec3vf<K> v2 = broadcast<vfloat<K>>(tri.v2,i);
          pre.intersectK(valid_i,ray,v0,v1,v2,UVIdentity<K>(),IntersectKEpilogM<M,K,filter>(ray,context,tri.geomID(),tri.primID(),i));
        }
      }

      /*! Test for K rays if they are occluded by any of the M triangles. */
      static __forceinline vbool<K> occluded(const vbool<K>& valid_i, Precalculations& pre, RayK<K>& ray, RayQueryContext* context, const Primitive& tri)
      {
        vbool<K> valid0 = valid_i;

        for (size_t i=0; i<M; i++)
        {
          if (!tri.valid(i)) break;
          STAT3(shadow.trav_prims,1,popcnt(valid_i),K);
          const Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0,i);
          const Vec3vf<K> v1 = broadcast<vfloat<K>>(tri.v1,i);
          const Vec3vf<K> v2 = broadcast<vfloat<K>>(tri.v2,i);
          pre.intersectK(valid0,ray,v0,v1,v2,UVIdentity<K>(),OccludedKEpilogM<M,K,filter>(valid0,ray,context,tri.geomID(),tri.primID(),i));
          if (none(valid0)) break;
        }
        return !valid0;
      }
      
      /*! Intersect a ray with M triangles and updates the hit. */
      static __forceinline void intersect(Precalculations& pre, RayHitK<K>& ray, size_t k, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(normal.trav_prims,1,1,1);
        pre.intersect(ray,k,tri.v0,tri.v1,tri.v2,UVIdentity<M>(),Intersect1KEpilogM<M,K,filter>(ray,k,context,tri.geomID(),tri.primID()));
      }

      /*! Test if the ray is occluded by one of the M triangles. */
      static __forceinline bool occluded(Precalculations& pre, RayK<K>& ray, size_t k, RayQueryContext* context, const Primitive& tri)
      {
        STAT3(shadow.trav_prims,1,1,1);
        return pre.intersect(ray,k,tri.v0,tri.v1,tri.v2,UVIdentity<M>(),Occluded1KEpilogM<M,K,filter>(ray,k,context,tri.geomID(),tri.primID()));
      }
    };
  }
}
