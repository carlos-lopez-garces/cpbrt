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