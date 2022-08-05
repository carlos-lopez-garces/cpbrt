#include "bssrdf.h"
#include "interpolation.h"
#include "parallel.h"
#include "scene.h"

// Computes the 1st moment of the Fresnel reflectance function Fr. The ith Fresnel moment
// is an integral involving the Fresnel reflectance function Fr over the hemisphere of 
// directions. It is approximated here using a polynomial of degree 5 (a Taylor polynomial?).
Float FresnelMoment1(Float eta) {
    Float eta2 = eta * eta;
    Float eta3 = eta2 * eta;
    Float eta4 = eta3 * eta;
    Float eta5 = eta4 * eta;

    if (eta < 1) {
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 + 2.49277f * eta4 - 0.68441f * eta5;
    }
    return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 - 1.27198f * eta4 + 0.12746f * eta5;
}

// Computes the 2nd moment of the Fresnel reflectance function Fr.
Float FresnelMoment2(Float eta) {
    Float eta2 = eta * eta;
    Float eta3 = eta2 * eta;
    Float eta4 = eta3 * eta;
    Float eta5 = eta4 * eta;

    if (eta < 1) {
        return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 + 0.07883f * eta4 + 0.04860f * eta5;
    }
    
    Float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
    return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 + 458.843f * r_eta + 404.557f * eta - 189.519f * eta2 + 54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
}

Float BeamDiffusionMS(Float sigma_s, Float sigma_a, Float g, Float eta, Float r) {
    const int nSamples = 100;
    Float Ed = 0;
    Float sigmap_s = sigma_s * (1 - g);
    Float sigmap_t = sigma_a + sigmap_s;
    Float rhop = sigmap_s / sigmap_t;
    Float D_g = (2 * sigma_a + sigmap_s) / (3 * sigmap_t * sigmap_t);
    Float sigma_tr = std::sqrt(sigma_a / D_g);
    Float fm1 = FresnelMoment1(eta), fm2 = FresnelMoment2(eta);
    Float ze = -2 * D_g * (1 + 3 * fm2) / (1 - 2 * fm1);

    Float cPhi = .25f * (1 - 2 * fm1), cE = .5f * (1 - 3 * fm2);
    for (int i = 0; i < nSamples; ++i) {
        Float zr = -std::log(1 - (i + .5f) / nSamples) / sigmap_t;
        Float zv = -zr + 2 * ze;
        Float dr = std::sqrt(r * r + zr * zr), dv = std::sqrt(r * r + zv * zv);
        Float phiD = Inv4Pi / D_g * (std::exp(-sigma_tr * dr) / dr - std::exp(-sigma_tr * dv) / dv);

        Float EDn = Inv4Pi * (zr * (1 + sigma_tr * dr) * std::exp(-sigma_tr * dr) / (dr * dr * dr) - zv * (1 + sigma_tr * dv) * std::exp(-sigma_tr * dv) / (dv * dv * dv));

        Float E = phiD * cPhi + EDn * cE;
        Float kappa = 1 - std::exp(-2 * sigmap_t * (dr + zr));
        Ed += kappa * rhop * rhop * E;
    }
    return Ed / nSamples;
}

Float BeamDiffusionSS(Float sigma_s, Float sigma_a, Float g, Float eta, Float r) {
    Float sigma_t = sigma_a + sigma_s, rho = sigma_s / sigma_t;
    Float tCrit = r * std::sqrt(eta * eta - 1);
    Float Ess = 0;
    const int nSamples = 100;
    for (int i = 0; i < nSamples; ++i) {
        Float ti = tCrit - std::log(1 - (i + .5f) / nSamples) / sigma_t;

        Float d = std::sqrt(r * r + ti * ti);
        Float cosThetaO = ti / d;

        Ess += rho * std::exp(-sigma_t * (d + tCrit)) / (d * d) * PhaseHG(cosThetaO, g) * (1 - FrDielectric(-cosThetaO, 1, eta)) * std::abs(cosThetaO);
    }
    return Ess / nSamples;
}

void ComputeBeamDiffusionBSSRDF(Float g, Float eta, BSSRDFTable *t) {
    t->opticalRadiusSamples[0] = 0;
    t->opticalRadiusSamples[1] = 2.5e-3f;
    for (int i = 2; i < t->mOpticalRadiusSamples; ++i) {
        t->opticalRadiusSamples[i] = t->opticalRadiusSamples[i - 1] * 1.2f;
    }

    for (int i = 0; i < t->nRhoSamples; ++i) {
        t->rhoSamples[i] = (1 - std::exp(-8 * i / (Float)(t->nRhoSamples - 1))) / (1 - std::exp(-8));
    }
        
    ParallelFor([&](int i) {
        for (int j = 0; j < t->mOpticalRadiusSamples; ++j) {
            Float rho = t->rhoSamples[i], r = t->opticalRadiusSamples[j];
            t->profile[i * t->mOpticalRadiusSamples + j] = 2 * Pi * r * (BeamDiffusionSS(rho, 1 - rho, g, eta, r) + BeamDiffusionMS(rho, 1 - rho, g, eta, r));
        }

        t->effectiveRho[i] = IntegrateCatmullRom(
            t->mOpticalRadiusSamples,
            t->opticalRadiusSamples.get(),
            &t->profile[i * t->mOpticalRadiusSamples],
            &t->profileCDF[i * t->mOpticalRadiusSamples]
        );
    }, t->nRhoSamples);
}

Spectrum SeparableBSSRDF::Sample_S(
    const Scene &scene, Float u1, const Point2f &u2, MemoryArena &arena, SurfaceInteraction *pi, Float *pdf
) const {
    // The spatial component can be sampled separately from the directional component. A point
    // on the surface/boundary is also sampled and returned in pi.
    Spectrum Sp = Sample_Sp(scene, u1, u2, arena, pi, pdf);
    if (!Sp.IsBlack()) {
        // The BxDF of the material at the sampled point pi is used to sample the directional component of S.
        // The SeparableBxDF may be an arbitrary BxDF.
        pi->bsdf = ARENA_ALLOC(arena, BSDF)(*pi);
        pi->bsdf->Add(ARENA_ALLOC(arena, SeparableBxDF)(this));
        pi->wo = Vector3f(pi->shading.n);
    }

    return Sp;
}

Spectrum SeparableBSSRDF::Sample_Sp(
    const Scene &scene, Float u1, const Point2f &u2, MemoryArena &arena, SurfaceInteraction *pi, Float *pdf
) const {
    // Choose projection axis for probe ray. The projection axis may be any of the basis vectors of
    // shading space (one of which, the normal n_o, is perpendicular to the assumed planar local surface,
    // and the other 2 are tangential to p_0 and parrallel to the assumed planar local surface). Since
    // the perpendicular projection probe ray is sure to intersect the surface, we choose it with 0.5
    // probability, whereas each of the 2 parallel projection probe rays get 0.25 probability each (they
    // are more likely to miss the surface, especially if the local surface is indeed planar or
    // convex).
    //
    // vz will be the chosen axis.
    Vector3f vx, vy, vz;
    if (u1 < 0.5f) {
        // Choose the perpendicular projection probe ray with 0.5 probability.
        vx = ss;
        vy = ts;
        vz = Vector3f(ns);
        u1 *= 2;
    } else if (u1 < 0.75f) {
        // Choose one of the parallel projection probe ray with 0.25 probability.
        vx = ts;
        vy = Vector3f(ns);
        vz = ss;
        u1 = (u1 - .5f) * 4;
    } else {
        // Choose the other parallel projection probe ray.
        vx = Vector3f(ns);
        vy = ss;
        vz = ts;
        u1 = (u1 - .75f) * 4;
    }

    // Choose spectral channel. Since the mean free path (the average distance that light travels through
    // this scattering medium before hitting a particle, i.e. before a scattering event occurs) varies
    // spectrally (i.e. across wavelengths), we choose at random a spectral channel; other Sample_Sp calls
    // will sample other channels (at random).
    int ch = Clamp((int) (u1 * Spectrum::nSamples), 0, Spectrum::nSamples - 1);
    u1 = u1 * Spectrum::nSamples - ch;

    // Sample a polar coordinate (r, phi) to sample, in turn, the radial scattering profile Sr.
    Float r = Sample_Sr(ch, u2[0]);
    if (r < 0) {
        // No scattering from this channel.
        return Spectrum(0.f);
    }
    Float phi = 2 * Pi * u2[1];

    // Cap the length of the sampling radius. This maximum radius length is also obtained from the radial
    // scattering profile Sr and its value is such that the corresponding sphere contains all the possible
    // points of incidence pi that are within reach of po before energy falls off completely. In other words,
    // incident radiance that enters the surface through points pi that are beyond this radius won't come out
    // through po.
    Float rMax = Sample_Sr(ch, 0.999f);
    if (r > rMax) {
        // The sampled radius places the point of incidence pi beyond reach. Radiance coming in through such
        // a point won't come out through po.
        return Spectrum(0.f);
    }

    // Compute probe ray. The point of intersection of this ray and the surface will be the sample incident
    // point pi. The origin of this ray lies on the surface of a sphere around po of radius rMax. The ray
    // has length probeRayLength and exits the sphere at point pTarget. Its direction is the chosen axis
    // of projection, vz.
    //
    // probeRayLength follows from the Pythagorean theorem and corresponds to the length of a chord of the
    // sphere of radius rMax. The chord joins the origin of the probe ray and the point where the probe
    // ray exits the sphere. Half of the chord is the opposite leg of the right trangle whose hypotenuse is
    // rMax and its adjacent leg is the sampled radius r.
    Float probeRayLength = 2 * std::sqrt(rMax * rMax - r * r);
    Interaction base;
    // The sampled polar coordinate (r, phi) and the chosen axis of projection vz are used to construct
    // the probe ray. The probe ray origin starts at po.
    base.p = po.p;
    // Then the origin moves to the sampled polar coordinate point.
    Vector3f sampledPolarDirection = r*(vx*std::cos(phi) + vy*std::sin(phi));
    base.p += sampledPolarDirection;
    // Then the origin moves along the direction of the chosen axis of projection, placing it on the
    // boundary of the sphere of radius rMax.
    base.p -= (probeRayLength*0.5) * vz;
    base.time = po.time;
    // The probe ray points from a point on the surface of the sphere of radius rMax to another point on
    // the surface of the sphere in the direction of the chosen axis of projection. Recall that probeRayLength
    // is the length of the chord connecting these 2 points. 
    Point3f pTarget = base.p + (probeRayLength * vz);

    // Intersect probe ray against the scene geometry within the sphere of radius rMax but along the chord.
    // There might be multiple intersections along the probe ray, which we store in an IntersectionChain (a linked list).
    struct IntersectionChain {
        SurfaceInteraction si;
        IntersectionChain *next = nullptr; 
    };
    IntersectionChain *chain = ARENA_ALLOC(arena, IntersectionChain)();

    // Accumulate chain of intersections along probe ray.
    IntersectionChain *ptr = chain;
    int nFound = 0;
    while (scene.Intersect(base.SpawnRayTo(pTarget), &ptr->si)) {
        // ptr->si is the intersection.
        base = ptr->si;
        // Is the intersection point on this surface or on other object's surface? Collect only the intersections
        // on this surface.
        if (ptr->si.primitive->GetMaterial() == material) {
            // The intersection is indeed on this surface. Store it in the list and keep going.
            IntersectionChain *next = ARENA_ALLOC(arena, IntersectionChain)();
            ptr->next = next;
            ptr = next;
            nFound++;
        }
    }

    // At this point we may have multiple points of incidence pi. Randomly choose one of the several intersections
    // as the sampled point of incidence pi.
    if (nFound == 0) {
        return Spectrum(0.0f);
    }
    int selected = Clamp((int)(u1 * nFound), 0, nFound - 1);
    while (selected-- > 0) {
        chain = chain->next;
    }
    *pi = chain->si;

    // Compute sample PDF and return the spatial BSSRDF term Sp.
    *pdf = Pdf_Sp(*pi) / nFound;
    return Sp(*pi);
}

Float SeparableBSSRDF::Pdf_Sp(const SurfaceInteraction &pi) const {
    // Express po-pi and ni (normal at pi) with respect to local coordinates at po.
    Vector3f piToPo = po.p - pi.p;
    Vector3f piToPoLocal(Dot(ss, piToPo), Dot(ts, piToPo), Dot(ns, piToPo));
    Normal3f niLocal(Dot(ss, pi.n), Dot(ts, pi.n), Dot(ns, pi.n));

    // Compute BSSRDF profile radius under projection along each axis.
    Float rProj[3] = { 
        std::sqrt(piToPoLocal.y * piToPoLocal.y + piToPoLocal.z * piToPoLocal.z),
        std::sqrt(piToPoLocal.z * piToPoLocal.z + piToPoLocal.x * piToPoLocal.x),
        std::sqrt(piToPoLocal.x * piToPoLocal.x + piToPoLocal.y * piToPoLocal.y)
    };

    // Return combined probability from all BSSRDF sampling strategies.
    Float pdf = 0, axisProb[3] = { .25f, .25f, .5f };
    Float chProb = 1/ (Float)Spectrum::nSamples;
    for (int axis = 0; axis < 3; ++axis) {
        for (int ch = 0; ch < Spectrum::nSamples; ++ch) {
            pdf += Pdf_Sr(ch, rProj[axis]) * std::abs(niLocal[axis]) * chProb * axisProb[axis];
        }
    }
    
    return pdf;
}

// Computes the radial profile of BSSRDF using spline-based interpolation of the tabulated
// samples of r and rho.
Spectrum TabulatedBSSRDF::Sr(Float r) const {
    Spectrum Sr(0.f);

    // Do spline-based interpolation for each of the channels of the spectrum.
    for (int ch = 0; ch < Spectrum::nSamples; ++ch) {
        // Reduce the dimensionality of Sr(eta, g, rho, sigma_t, r) by fixing sigma_t=1 and
        // bundling sigma_t and r into a unitless optical radius r_optical = sigma_t * r. 
        // This is effectively a change of variable that requires a scaling factor (a Jacobian
        // determinant): Sr(eta, g, rho, sigma_t, r) = sigma_t^2 * Sr(eta, g, rho, 1, r_optical).
        Float rOptical = r * sigma_t[ch];

        // Compute spline weights to interpolate BSSRDF on channel ch. The weights are the
        // coefficients of a polynomial interpolator of a function f (the Sr table in this case):
        //
        // p(x) = w0*f(x_-1) + w1*f(x_0) + w2*f(x_1) + w3*f(x_2).
        //
        // An interpolator is obtained for rho and another for optical radius.

        // An offset is the index (into the row) that corresponds to the start endpoint of the
        // subinterval (a pair of samples or r or rho) that contains the lookup value (input r
        // and TabulatedBSSRDF::rho).  
        int rhoOffset;
        int radiusOffset;
        // Interpolator weights for rho and optical radius.
        Float rhoWeights[4];
        Float radiusWeights[4];
        if (
            !CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch], &rhoOffset, rhoWeights)
            || !CatmullRomWeights(table.mOpticalRadiusSamples, table.opticalRadiusSamples.get(), rOptical, &radiusOffset, radiusWeights)
        ) {
            // rho or r are not in the domain of the tabulated radial profile Sr.
            continue;
        }

        // Set BSSRDF value Sr[ch] using tensor spline interpolation.
        Float sr = 0;
        // 4 weights.
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Float weight = rhoWeights[i] * radiusWeights[j];
                if (weight != 0) {
                    sr += weight * table.EvalProfile(rhoOffset + i, radiusOffset + j);
                }
            }
        }

        // Cancel marginal PDF factor of 2*Pi*r_optical from tabulatd BSSRDF profile. It
        // was introduced when evaluating the tabulated profile and is there for facilitating
        // importance sampling, but it's not part of the definition of Sr.
        if (rOptical != 0) {
            sr /= 2 * Pi * rOptical;
        }

        Sr[ch] = sr;
    }

    // Transform unitless value into world space units. The Sr value thus far computed is unitless
    // as a result of a change of variable r_optical = sigma_t * r that reduces the dimensionality
    // of Sr (see above).
    Sr *= sigma_t * sigma_t;

    return Sr.Clamp();
}

Float TabulatedBSSRDF::Sample_Sr(int ch, Float u) const {
    if (sigma_t[ch] == 0) {
        // Radiance is instantly attenuated or extinguished at wavelength ch. No radius to return.
        return -1;
    }

    // TODO: implement SampleCatmullRom2D. 
    Float rOptical = SampleCatmullRom2D(
        table.nRhoSamples,
        table.mOpticalRadiusSamples,
        table.rhoSamples.get(),
        table.opticalRadiusSamples.get(),
        table.profile.get(),
        table.profileCDF.get(),
        rho[ch],
        u
    );

    // Radius is optical radius divided by extinction coefficient, r = r_optical / sigma_t.
    return rOptical / sigma_t[ch];
}

Float TabulatedBSSRDF::Pdf_Sr(int ch, Float r) const {
    Float rOptical = r * sigma_t[ch];

    int rhoOffset, radiusOffset;
    Float rhoWeights[4], radiusWeights[4];
    if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch], &rhoOffset, rhoWeights)
        || !CatmullRomWeights(table.mOpticalRadiusSamples, table.opticalRadiusSamples.get(), rOptical, &radiusOffset, radiusWeights)) {
        return 0.f;
    }

    Float sr = 0, rhoEff = 0;
    for (int i = 0; i < 4; ++i) {
        if (rhoWeights[i] == 0) {
            continue;
        }

        rhoEff += table.effectiveRho[rhoOffset + i] * rhoWeights[i];
        for (int j = 0; j < 4; ++j) {
            if (radiusWeights[j] == 0) {
                continue;
            }

            sr += table.EvalProfile(rhoOffset + i, radiusOffset + j) * rhoWeights[i] * radiusWeights[j];
        }
    }

    if (rOptical != 0) {
        sr /= 2 * Pi * rOptical;
    }
    return std::max((Float)0, sr * sigma_t[ch] * sigma_t[ch] / rhoEff);
}