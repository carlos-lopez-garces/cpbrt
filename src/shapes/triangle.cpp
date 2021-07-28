#include "triangle.h"

TriangleMesh::TriangleMesh(
    const Transform &ObjectToWorld,
    int nTriangles,
    const int *vertexIndices,
    int nVertices,
    const Point3f *P,
    const Vector3f *S,
    const Normal3f *N,
    const Point2f *UV,
    const std::shared_ptr<Texture<Float>> &alphaMask,
    const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
    const int *fIndices
) : nTriangles(nTriangles),
    nVertices(nVertices),
    vertexIndices(vertexIndices, vertexIndices + 3*nTriangles),
    alphaMask(alphaMask),
    shadowAlphaMask(shadowAlphaMask) {

    // Copy input arrays into class members, transforming their values to world space.
    // Ray intersection thus doesn't have to transform the ray to object space; the test
    // can be done in world space.

    p.reset(new Point3f[nVertices]);
    for (int i = 0; i < nVertices; ++i) {
        p[i] = ObjectToWorld(P[i]);
    }

    if (UV) {
        uv.reset(new Point2f[nVertices]);
        memcpy(uv.get(), UV, nVertices * sizeof(Point2f));
    }

    if (N) {
        n.reset(new Normal4f[nVertices]);
        for (int i = 0; i < nVertices; ++i) {
            // Transforms have an overloaded () operator for transforming normals.
            n[i] = ObjectToWorld(N[i]);
        }
    }

    // Per-vertex tangent vectors.
    if (S) {
        s.reset(new Vector3f[nVertices]);
        for (int i = 0; i < nVertices; ++i) {
            s[i] = ObjectToWorld(S[i]);
        }
    }

    if (fIndices) {
        faceIndices = std::vector<int>(fIndices, fIndices + nTriangles);
    }
}

Bounds3f Triangle::ObjectBound() const {
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    // Bounding box that encloses the 3 vertices, in object space.
    return Union(
        Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
        (*WorldToObject)(p2)
    );
}

Bounds3f Triangle::WorldBound() const {
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    // Bounding box that encloses the 3 vertices, in world space.
    return Union(Bounds3f(p0, p1), p2);
}

bool Triangle::Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *si, bool testAlphaTexture) const {
    // Get triangle vertices.
    const Point3f &p0 = mesh->p[v[0]];
    const Point3f &p1 = mesh->p[v[1]];
    const Point3f &p2 = mesh->p[v[2]];

    // Perform ray-triangle intersection test. This test is done in a "ray-triangle intersection
    // coordinate system". Vectors are transformed to this coordinate system with the affine matrix
    // transformation M=SPT, where T is a translation, P is a coordinate permutation, and S is a shear.

    // The translation places the origin of the coordinate system at the ray's origin. The matrix
    // multiplication by T is equivalent to a subtraction of the ray's origin from the point/vector.
    //
    // Point-point subtraction doesn't make sense, that's why the operands are a point and a vector
    // (the subtraction displaces the point in the direction of the vector, by the magnitude of the
    // vector). NOTE: These are not homogeneous coordinates.
    Point3f p0t = p0 - Vector3f(ray.o);
    Point3f p1t = p1 - Vector3f(ray.o);
    Point3f p2t = p2 - Vector3f(ray.o);

    // The coordinate permutation P rotates the components of the permuted vector in a way such that
    // their z component ends up with the value of the dimension that is largest in absolute value
    // in the ray direction vector.
    //
    // kx, ky, kz are the component numbers that will end up in x, y, and z, respectively, not the
    // actual coordinate values.
    int kz = MaxDimension(Abs(ray.d));
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    Vector3f d = Permute(ray.d, kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);

    // The shear S takes the ray direction vector to the k basis vector (0,0,1) of the coordinate
    // system. The +z axis of the coordinate system will point in that direction. Note that we defer
    // transforming the z coordinate of the vertices until later: if the ray doesn't intersect the
    // coordinate system's xy-plane at a point inside the triangle (the test of which needs only the
    // transformed x and y coordinates), there's no point in transforming the z coordinate now. 
    // TODO: compute the Sx, Sy, Sz constants in the Ray constructor, instead of recomputing them
    // on every intersection test.
    Float Sx = -d.x / d.z;
    Float Sy = -d.y / d.z;
    Float Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;

    // At this point, the x and y coordinates of the 3 vertices of the triangle have coordinates in
    // the intended coordinate system. The z coordinate is transformed only if the following test passes.

    // Compute edge function coefficients e0, e1, and e2. The signed edge function e(p) tells where
    // the point p lies with respect to the line where the given edge lies: to the right, left, or
    // on it (you have to view the edge vector from tail to head for this relative orientation to
    // make sense). The signed edge function e(p) computes twice the area of the triangle formed by a
    // given edge of our triangle and the point we are testing. Now, consider points P and Q, each on
    // either side of the edge, but at the same distance from it: the absolute value of the double area
    // will be the same for P and Q (|e(P)| = |e(Q)|), but the signs will be opposite.
    //
    // Furthermore, if the edges of the triangle are sequenced like this p0p1, p1p2, p2p0, then the point
    // p will be on the same side of each of the 3 edges when p is inside the triangle (if it's outside
    // 2 of the edges will have one sign, and the other the opposite).
    //
    // The expressions of the 3 signed edge functions are derived from that of the magnitude of the
    // cross product between the edge's vector and the tail of the edge vector and p, which gives the
    // area of the parallelogram that these vectors define, which is exactly twice the area of the triangle.
    // There aren't any squares like you'd imagine, because the cross product of 2 vectors that lie on the
    // xy-plane (like the vectors of our points here) has 0 x and y components: the square root of the square
    // of the z coordinate is the z coordinate (Euclidian distance formula).
    //
    // The point p we are testing here is always (0,0) because that's the projection of the ray onto the
    // xy-plane of the ray-triangle intersection coordinate system. The (x,y) coordinates of p0t, p1t, and
    // p2t (the vertices of the triangle transformed to this coordinate system) are the projections of
    // these points onto the xy-plane and the vectors between them define the edges of each of the signed
    // edge functions.
    Float e0 = p1t.x*p2t.y - p1t.y*p2t.x;
    Float e1 = p2t.x*p0t.y - p2t.y*p0t.x;
    Float e2 = p0t.x*p1t.y - p0t.y*p1t.x;
    if (sizeof(Float) == sizeof(float) && (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        // The ray hit one of the edges of the triangle exactly or a shared vertex. But increasing the
        // floating-point precision might reveal that that wasn't the case, actually, and that it was the
        // result of rounding a non-representable real-number result up or down to 0. Reevaluate the signed
        // edge functions with double precision.
        double p2txp1ty = (double) p2t.x * (double) p1t.y;
        double p2typ1tx = (double) p2t.y * (double) p1t.x;
        e0 = (float) (p2typ1tx - p2txp1ty);
        double p0txp2ty = (double) p0t.x * (double) p2t.y;
        double p0typ2tx = (double) p0t.y * (double) p2t.x;
        e1 = (float) (p0typ2tx - p0txp2ty);
        double p1txp0ty = (double) p1t.x * (double) p0t.y;
        double p1typ0tx = (double) p1t.y * (double) p0t.x;
        e2 = (float) (p1typ0tx - p1txp0ty);
    }
    // For the point to be inside the triangle (and the ray to hit it), it has to lie on the same side of
    // all edges: either to the left of all of them or to the right.
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0)) {
        // The signed edge functions differ in sign: the point is outside the triangle for sure.
        return false;
    }
    // TODO: explain.
    Float determinant = e0 + e1 + e2;
    if (determinant == 0) {
        // The triangle is orthogonal to the coordinate system's xy-plane, so there's no way the ray can
        // hit it.
        return false;
    }

    // At this point we know that the ray is aligned with the interior of the triangle, but we don't know
    // if it is within the parametric range of the ray. And we haven't computed the actual intersection
    // point either; neither in the ray-triangle coordinate system nor in world space.

    // Apply the shear transformation to the z coordinate (we had deferred it).
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;

    // Compute scaled hit distance to triangle and test against ray t range.
    // TODO: explain.
    Float tScaled = e0*p0t.z + e1*p1t.z + e2*p2t.z;
    if (determinant < 0 && (tScaled >= 0 || tScaled < ray.tMax * determinant)) {
        return false;
    } else if (determinant > 0 && (tScaled <= 0 || tScaled > ray.tMax * determinant)) {
        return false;
    }

    // The barycentric coordinates (b0, b1, b2) of the intersection point are used to compute the
    // coordinates of the intersection point. A linear combination of the world-space coordinates of the
    // vertices p0, p1, p2 of the triangle interpolates them to give the world-space coordinates of the
    // intersection point; the coefficients of the linear combination are the barycentric coordinates:
    //
    // pHit = b0p0  + b1p1 + b2p2
    //
    // The barycentric coordinates result from normalizing the values of the signed edge functions:
    //
    // (b1, b2, b3) = (e0/(e0+e1+e2), e1/(e0+e1+e2), e2/(e0+e1+e2))
    //
    // So b1 + b2 + b3 = 1.
    //
    // Recall that signed edge function values are the signed areas of the 3 parallelograms defined by
    // the 3 edges and the intersection point. The areas of the parallelograms are twice the areas of
    // the triangles defined by the 3 edges and the intersection point, which means that these values
    // are proportional to the areas of the triangles. The normalization of the parallelogram areas
    // and the normalization of the triangle areas result in exactly the same coordinates:
    //
    // (b1, b2, b3)
    // = (e0/(e0+e1+e2), e1/(e0+e1+e2), e2/(e0+e1+e2)) 
    // = (2et0/2(et0+et1+et2), 2et1/2(et0+et1+et2), 2et2/2(et0+et1+et2))
    // = (et0/(et0+et1+et2), et1/(et0+et1+et2), et2/(et0+et1+et2))        (Factor 2 cancels out.)
    //
    // which is exactly the general definition of barycentric coordinates (also known as areal
    // coordinates).
    //
    // (et denotes the area of 1 of the component triangles of the corresponding parallelogram.)
    //
    // The sum e0+e1+e2 is already held by the determinant variable.
    Float invDeterminant = 1 / determinant;
    Float b0 = e0 * invDeterminant;
    Float b1 = e1 * invDeterminant;
    Float b2 = e2 * invDeterminant;

    // Compute ray t value for intersection.
    Float t = tScaled * invDeterminant;

    // Ensure that computed triangle t is conservatively greater than zero.
    // Compute deltaZ term for triangle t error bounds.
    Float maxZt = MaxComponent(Abs(Vector3f(p0t.z, p1t.z, p2t.z)));
    Float deltaZ = gamma(3) * maxZt;
    // Compute deltaX and deltaY terms for triangle t error bounds.
    Float maxXt = MaxComponent(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
    Float maxYt = MaxComponent(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);
    // Compute deltaE term for triangle t error bounds.
    Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
    // Compute deltaT term for triangle t error bounds and check t.
    Float maxE = MaxComponent(Abs(Vector3f(e0, e1, e2)));
    Float deltaT = 3 * (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * std::abs(invDeterminant);
    if (t <= deltaT) {
        return false;
    }

    // Compute triangle partial derivatives.
    Vector3f dpdu;
    Vector3f dpdv;
    // (u,v) parameterization of the 3 triangle vertices.
    Point2f uv[3];
    GetUVs(uv);
    // From vertex p2, compute the vectors to the other 2 vertices in parametric (u,v) space and
    // world space; referred to as deltas (differences).
    // Compute (u,v) deltas.
    Vector2f duv02 = uv[0] - uv[2];
    Vector2f duv12 = uv[1] - uv[2];
    // Compute world space p deltas.
    Vector3f dp02 = p0 - p2;
    Vector3f dp12 = p1 - p2;
    // Express the partial derivatives of p with respect to u and v with a matrix equation. This
    // matrix equation comes from expressing the coordinates of a vertex of the triangle in terms
    // of the (u,v) parameters and the partial derivatives ...
    //
    // [u0-u2  v0-v2] [partialP/partialU] = [p0-p2]
    // [u1-u2  v1-v2] [partialP/partialV] = [p1-p2]
    //
    //... and then solving for the partial derivatives. Solving for the partial derivatives involves 
    // inverting the square matrix of (u,v) deltas:
    //
    //                                     -1
    // [partialP/partialU] = [u0-u2  v0-v2]  [p0-p2]
    // [partialP/partialV]   [u1-u2  v1-v2]  [p1-p2]
    //
    // You don't need to do row reduction of the (u,v) deltas matrix augmented with the identity
    // matrix to obtain the inverse: square matrices can easily be inverted with a closed-form
    // expression based on the determinant.
    Float deltasDeterminant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    if (deltasDeterminant == 0) {
        // The (u,v) deltas matrix is not invertible. Obtain an arbitrary orthonormal set of vectors
        // about the triangle's normal for the partial derivatives.
        // TODO: why?
        CoordinateSystem(Normalize(Cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);
    } else {
        // Invert the (u,v) deltas matrix and solve for the partial derivatives matrix.
        Float invDeltasDeterminant = 1 / deltasDeterminant;
        dpdu = ( duv12[1] * dp02 - duv02[1] * dp12) * invDeltasDeterminant;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invDeltasDeterminant;
    }

    // Compute error bounds for intersection.
    Float xAbsSum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
    Float yAbsSum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
    Float zAbsSum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
    Vector3f pError = gamma(7) * Vector3f(xAbsSum, yAbsSum, zAbsSum);

    // Interpolate the (u,v) parametric coordinates and world-space coordinates of the triangle
    // vertices to obtain (u,v) and world-space coordinates for the intersection point. Use the
    // barycentric coordinates of the intersection point to do the interpolation.
    Point3f pHit = b0*p0 + b1*p1 + b2*p2;
    Point2f uvHit = b0*uv[0] + b1*uv[1] + b2*uv[2];

    // TODO: Test intersection against alpha texture, if present.

    // Fill in SurfaceInteraction. The pair of Normal3f(0,0,0) inputs are dndu and dndv: the partial
    // derivatives of the normal with respect to (u,v). Since the triangle is flat, the normal doesn't
    // change across its surface as (u,v) varies over the domain.
    *si = SurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv, Normal3f(0,0,0), Normal3f(0,0,0), ray.time, this);
    // Compute geometric normal, using triangle winding order consistently across invocations and across
    // triangles of the mesh.
    si->n = si->shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
    // Ensure correct orientation of the geometric and shading normals.
    if (mesh->n) {
        // The mesh uses per-vertex normals.
        si->n = Faceforward(si->n, si->shading.n);
    } else if (reverseOrientation ^ transformSwapsHandedness) {
        // Reverse both. (The condition is either-or.)
        si->n = si->shading.n = -si->n;
    }

    if (mesh->n || mesh->s) {
        // The mesh uses per-vertex normals or tangents. Initialize shading geometry.

        Normal3f ns;
        if (mesh->n) {
            // Compute shading normal by interpolating the mesh's per-vertex normals using the barycentric
            // coordinates of the intersection point.
            ns = Normalize(b0*mesh->n[v[0]] + b1*mesh->n[v[1]] + b2*mesh->n[v[2]]);
        } else {
            // The mesh doesn't include per-vertex normals. Use geometric normal for shading.
            ns = si->n;
        }

        Vector3f ss;
        if (mesh->s) {
            // Compute shading tangent by interpolating the mesh's per-vertex tangents using the barycentric
            // coordinates of the intersection point.
            ss = Normalize(b0*mesh->s[v[0]] + b1*mesh->s[v[1]] + b2*mesh->s[v[2]]);
        } else {
            // The mesh doesn't include per-vertex tangents. Use one of the partial derivatives; they are
            // tangent to the surface and orthogonal to the normal.
            ss = Normalize(si->dpdu);
        }

        // Compute shading bitangent. The bitangent vector completes the orthonormal frame formed by the
        // normal and tangent vectors. It's also tangent to the surface.
        Vector3f ts = Cross(ns, ss);
        if (ts.LengthSquared() > 0.f) {
            ts = Normalize(ts);
            // When the mesh includes per-vertex normals, the shading normal was the result of interpolating
            // these per-vertex normals. Recompute the shading tangent now that we have the other 2 vectors
            // of the orthonormal frame.
            ss = Cross(ts, ns);
        } else {
            // The shading normal is equal to the geometric normal, which accurately describes the orientation
            // of the triangle. Recompute the shading tangent and compute the shading bitangent using the shading
            // (geometric) normal.
            CoordinateSystem((Vector3f) ns, &ss, &ts);
        }
        if (reverseOrientation) ts = -ts;

        // Compute partial derivatives of the shading normal with respect to (u,v).
        Normal3f dndu, dndv;
        if (mesh->n) {
            Vector2f duv02 = uv[0] - uv[2];
            Vector2f duv12 = uv[1] - uv[2];
            Normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
            Normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
            Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            bool degenerateUV = std::abs(determinant) < 1e-8;
            if (degenerateUV) {
                Vector3f dn = Cross(Vector3f(mesh->n[v[2]] - mesh->n[v[0]]), Vector3f(mesh->n[v[1]] - mesh->n[v[0]]));
                if (dn.LengthSquared() == 0) {
                    dndu = dndv = Normal3f(0, 0, 0);
                }
                else {
                    Vector3f dnu, dnv;
                    CoordinateSystem(dn, &dnu, &dnv);
                    dndu = Normal3f(dnu);
                    dndv = Normal3f(dnv);
                }
            } else {
                Float invDet = 1 / determinant;
                dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
                dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
            }
        } else {
            dndu = dndv = Normal3f(0, 0, 0);
        }

        si->SetShadingGeometry(ss, ts, dndu, dndv, true);
    }

    *tHit = t;
    return true;
}

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    const Transform *ObjectToWorld, 
    const Transform *WorldToObject,
    bool reverseOrientation,
    int nTriangles,
    const int *vertexIndices,
    int nVertices,
    const Point3f *p,
    const Vector3f *s,
    const Normal3f *n,
    const Point2f *uv,
    const std::shared_ptr<Texture<Float>> &alphaMask,
    const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
    const int *faceIndices
) {
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
        *ObjectToWorld,
        nTriangles,
        vertexIndices,
        nVertices,
        p, s, n, uv,
        alphaMask,
        shadowAlphaMask,
        faceIndices
    );

    // The TriangleMesh constructor doesn't create the Triangle shapes.
    std::vector<std::shared_ptr<Shape>> triangles;
    triangles.reserve(nTriangles);
    for (int i = 0; i < nTriangles; ++i) {
        triangles.push_back(std::make_shared<Triangle>(
            // The Triangle shape uses the index i to find its data in the TriangleMesh.
            ObjectToWorld, WorldToObject, reverseOrientation, mesh, i
        ));
    }

    return triangles;
}