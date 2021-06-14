#include "scene.h"

bool Scene::Intersect(const Ray &ray, SurfaceInteraction *si) const {
    return aggregate->Intersect(ray, si);
}

bool Scene::IntersectP(const Ray &ray) const {
    return aggregate->IntersectP(ray);
}

 bool Scene::IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *si, Spectrum *transmittance) const {
     // TODO: implement.
     return false;
 }