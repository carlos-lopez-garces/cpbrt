#include "core/medium.h"
#include "core/memory.h"
#include "core/spectrum.h"
#include "core/transform.h"

class GridDensityMedium : public Medium {
private:
    // Absorption cross section. The probability that radiance is absorbed per unit
    // distance. Absorption is a source of attenuation.
    const Spectrum sigma_a;

    // Absorption cross section. The probability that radiance is absorbed per unit
    // distance. Absorption is a source of attenuation.
    const Spectrum sigma_s;

    // Attenuation or extinction coefficient. The combined effect of absorption and
    // outscattering: sigma_t(p,w) = sigma_a(p,w) + sigma_s(p,w).
    Float sigma_t;

    // Henyey-Greenstein phase function's asymmetry parameter. Negative g values describe
    // phase functions that primarily scatter light back in the incident direction, and
    // positive g values, forward in the direction that light is travelling. Isotropic
    // phase functions, which describe equal scattering in all directions, have g = 0.
    const Float g;

    // Dimensions of the grid.
    const int nx, ny, nz;

    // Density samples at the points of the grid. Used to reconstruct the volume density
    // function via interpolation. Samples may be obtained via physical simulation, CT
    // scan, etc.
    std::unique_ptr<Float[]> density;
    
    Float reciprocalMaxDensity;
    const Transform WorldToMedium;


public:
    GridDensityMedium(
        const Spectrum &sigma_a,
        const Spectrum &sigma_s,
        Float g,
        int nx, int ny, int nz,
        const Transform &mediumToWorld,
        const Float *d
    ) : sigma_a(sigma_a), 
        sigma_s(sigma_s), 
        g(g),
        nx(nx), ny(ny), nz(nz),
        WorldToMedium(Inverse(mediumToWorld)),
        density(new Float[nx*ny*nz]) {

        // Precompute values for Monte Carlo sampling.
        sigma_t = (sigma_a + sigma_s)[0];
        Float maxDensity = 0;
        for (int i = 0; i < nx*ny*nz; ++i) {
            maxDensity = std::max(maxDensity, density[i]);
        }
        reciprocalMaxDensity = 1 / maxDensity;
    }

    // Reconstructs the volume density function from the grid samples at the
    // given point via interpolation.
    Float Density(const Point3f &p) const;

    // Computes the density at the given integer sample position. Sample points
    // outside the grid are given a density of 0.
    Float D(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) {
            return 0;
        }
        return density[(p.z * ny + p.y) * nx + p.x];
    }

    Spectrum Sample(
        const Ray &ray,
        Sampler &sampler,
        MemoryArena &arena,
        MediumInteraction *mi
    ) const;
};