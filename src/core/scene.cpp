#include "scene.h"

bool Scene::Intersect(const Ray &ray, SurfaceInteraction *si) const {
    return aggregate->Intersect(ray, si);
}

bool Scene::IntersectP(const Ray &ray) const {
    return aggregate->IntersectP(ray);
}