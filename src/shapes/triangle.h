
#ifndef CPBRT_SHAPES_TRIANGLE_H
#define CPBRT_SHAPES_TRIANGLE_H

#include "core/shape.h"

struct TriangleMesh {
    const int nTriangles;
    const int nVertices;
    std::vector<int> vertexIndices;
    // World-space coordinates of (unique) vertices.
    std::unique_ptr<Point3f[]> p;
    // World-space per-vertex normals.
    std::unique_ptr<Normal3f[]> n;
    // Wolrd-space per-vertex tangent vectors.
    std::unique_ptr<Vector3f[]> s;
    // Mesh parameterization.
    std::unique_ptr<Point2f[]> uv;
    std::shared_ptr<Texture<Float>> alphaMask;
    std::shared_ptr<Texture<Float>> shadowAlphaMask;
    // ?
    std::vector<int> faceIndices;

    // Doesn't create the Triangle shapes, just the mesh (their connectivity info).
    TriangleMesh(
        const Transform &ObjectToWorld,
        int nTriangles,
        // Indices into the P array of vertices. The indices of a given triangle
        // are stored contiguously. 
        const int *vertexIndices,
        int nVertices,
        // Object-space coordinates of (unique) vertices.
        const Point3f *P,
        // Object-space per-vertex tangent vectors.
        const Vector3f *S,
        // Object-space per-vertex normals.
        const Normal3f *N,
        // Per-vertex parameterization.
        const Point2f *UV,
        const std::shared_ptr<Texture<Float>> &alphaMask,
        const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
        // ?
        const int *faceIndices
    );
};

class Triangle : public Shape {
private:
    // The Triangle class doesn't store the coordinates of its vertices, just a pointer
    // to the start of its 3 indices in the mesh.
    std::shared_ptr<TriangleMesh> mesh;
    // A pointer to the start of its 3 indices in TriangleMesh.vertexIndices.
    const int *v;
    int faceIndex;

public:
    Triangle(
        const Transform *ObjectToWorld,
        const Transform *WorldToObject,
        bool reverseOrientation,
        const std::shared_ptr<TriangleMesh> *mesh,
        int index
    ) : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        mesh(mesh) {

        v = &mesh->vertexIndices[3 * index];
        faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[index] : 0;
    }

    // Bounds in object space.
    Bounds3f ObjectBound() const;

    // Bounds in world space.
    Bounds3f WorldBound() const;

    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *si, bool testAlphaTexture = true) const;

private:
    // Puts the UVs in the input array.
    void GetUVs(Point2f uv[3]) const {
        if (mesh->uv) {
            uv[0] = mesh->uv[v[0]];
            uv[1] = mesh->uv[v[1]];
            uv[2] = mesh->uv[v[2]];
        } else {
            // Default parameterization. See the constructor: caller-supplied UVs are optional.
            uv[0] = Point2f(0, 0);
            uv[1] = Point2f(1, 0);
            uv[2] = Point2f(1,1);
        }
    }
};

// Creates a TriangleMesh and the component Triangle shapes. Returns the array of
// Triangle shapes (the TriangleMesh is referenced by each Triangle via a shared_ptr).
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
    const int *faceIndices = nullptr
);

#endif // CPBRT_SHAPES_TRIANGLE_H