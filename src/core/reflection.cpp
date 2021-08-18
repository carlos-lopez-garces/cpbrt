#include "cpbrt.h"
#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "scene.h"
#include "interaction.h"


// Evaluates the Fresnel reflectance equation between 2 dielectric media, assuming that light is
// unpolarized. cosThetaI is the angle of incidence measured from the normal; etaI is the refraction
// index of the medium that light is traveling through before reaching the interface with the
// other medium, whose refraction index etaT is. (A refraction index is a property of the medium:
// the ratio of the speed of light in a vacuum to the speed of light through the medium. Refraction
// indices for dielectrics are assumed to be real numbers.)
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
    cosThetaI = Clamp(cosThetaI, -1, 1);

    // Make sure that etaI is really the index of the incident medium and etaT the index of the transmitted
    // medium. cosThetaI was measured between the surface's normal and the direction vector of incidence: if
    // it's negative, the 2 are in opposite hemispheres, and the ray is hitting the surface of the medium
    // from within that surface's medium: this medium's refraction index should really be etaI and not etaT.
    if (cosThetaI <= 0.f) {
        // Incidence vector is not inside the transmission medium, but the incident one.
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute cosThetaT using Snell's law: etaI*sinThetaI = etaT*sinThetaT.
    // sinThetaI and cosThetaT are computed using a trigonometric identity: sin(theta)^2 + cos(theta)^2 = 1.
    Float sinThetaI = std::sqrt(std::max((Float) 0, 1 - cosThetaI*cosThetaI));
    Float sinThetaT = etaI / etaT * sinThetaI;
    if (sinThetaT >= 1) {
        // Total internal reflection: light grazes the boundary of a medium with lower refraction index.
        // Snell's law doesn't have a solution, so refraction can't occur: light gets reflected back into
        // the incident medium.
        return 1;
    }
    Float cosThetaT = std::sqrt(std::max((Float) 0, 1 - sinThetaT*sinThetaT));

    // Evaluate Fresnel reflectance equation for the parallel polarized component of light.
    Float parallelR = ((etaT*cosThetaI) - (etaI*cosThetaT)) / ((etaT*cosThetaI)+(etaI*cosThetaT));

    // Evaluate Fresnel reflectance equation for the perpendicular polarized component of light.
    Float perpendicularR = ((etaI*cosThetaI) - (etaT*cosThetaT)) / ((etaI*cosThetaI) + (etaT*cosThetaT));

    // The reflectance of unpolarized light is the average of the parallel and perpendicular polarized reflectances.
    return (parallelR*parallelR + perpendicularR*perpendicularR) / 2;
}

// Evaluates the general Fresnel equation, but assuming that the transmission medium is a dielectric.
// Refraction indices are spectra because they are wavelength-dependent. cosThetaI is the angle of
// incidence measured from the normal; etaI and etaT are the indices of refraction of the incident
// and transmission media, respectively; and k is the imaginary part of the index of refraction of
// the incident medium (assumed to be a conductor) also known as absorption coefficient.
Spectrum FrConductor(Float cosThetaI, const Spectrum &etaI, const Spectrum &etaT, const Spectrum &k) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    // Relative refraction index.
    Spectrum eta = etaT / etaI;
    Spectrum etaK = k / etaI;

    Float cosThetaI2 = cosThetaI * cosThetaI;
    Float sinThetaI2 = 1. - cosThetaI2;
    Spectrum eta2 = eta * eta;
    Spectrum etaK2 = etaK * etaK;

    // Compute terms.
    Spectrum t0 = eta2 - etaK2 - sinThetaI2;
    Spectrum a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etaK2);
    Spectrum t1 = a2plusb2 + cosThetaI2;
    Spectrum a = Sqrt(0.5f * (a2plusb2 + t0));
    Spectrum t2 = (Float) 2 * cosThetaI * a;
    Spectrum t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    Spectrum t4 = t2 * sinThetaI2;

    Spectrum perpendicularR = (t1 - t2) / (t1 + t2);
    Spectrum parallelR = perpendicularR * (t3 - t4) / (t3 + t4);

    // The reflectance of unpolarized light is the average of the parallel and perpendicular polarized reflectances.
    return (perpendicularR + parallelR) / 2;
}

Spectrum FresnelDielectric::Evaluate(Float cosThetaI) const {
    return FresnelDielectric(cosThetaI, etaI, etaT);
}

Spectrum FresnelConductor::Evaluate(Float cosThetaI) const {
    return FrConductor(std::abs(cosThetaI), etaI, etaT, k);
}

Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
    Spectrum r(0.);
    for (int i = 0; i < nSamples; ++i) {
        Vector3f wi;
        Float pdf = 0;
        Spectrum f = Sample_f(w, &wi, u[i], &pdf);
        if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
    }
    return r / nSamples;
}

Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
    Spectrum r(0.f);
    for (int i = 0; i < nSamples; ++i) {
        Vector3f wo, wi;
        wo = UniformSampleHemisphere(u1[i]);
        Float pdfo = UniformHemispherePdf(), pdfi = 0;
        Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
        if (pdfi > 0) {
            r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
        }
    }
    return r / (Pi * nSamples);
}

Spectrum BxDF::Sample_f(
    // wo and wi are expressed with respect to the local coordinate system.
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &u,
    Float *pdf,
    BxDFType *sampledType
) const {
    // The incident direction wi corresponding to the outgoing direction wo may be any
    // one on the hemisphere centered at the point u and on the side of the surface where
    // wo exits (which isn't always the side of the surface's normal). This incident direction
    // wi is sampled from a cosine-weighted distribution, which assigns a higher probability
    // to directions at the top of the hemisphere than those at the bottom (why is this
    // desirable?).
    //
    // This base implementation of Sample_f is not for BTDFs, which should sample wi from the
    // opposite hemisphere, nor is it for BRDFs with delta distribution (which map an outgoing
    // direction wo to a single fixed incident direction wi).
    *wi = CosineSampleHemisphere(u);

    if (wo.z < 0) {
        // CosineSampleHemisphere samples the hemisphere on the normal's side. When the outgoing
        // direction wo lies on the hemisphere opposite to the normal, wi will have been sampled
        // from the wrong hemisphere. Bring it to wo's side.
        wi->z *= -1;
    }

    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

// Computes the PDF with which the base implementation of Sample_f samples the incident direction
// wi. This PDF is of a cosine-weighted distribution: p(w)=r*cos(theta)/pi, where r=1 is the radius
// of the unit hemisphere and theta is measured from the hemisphere's axis.
Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

Spectrum ScaledBxDF::f(const Vector3f &wo, const Vector3f &wi) const {
    // The product of 2 spectrums is sample-wise.
    return scale * bxdf->f(wo, wi);
}

Spectrum ScaledBxDF::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &sample,
    Float *pdf,
    BxDFType *sampledType
) const {
    // The product of 2 spectrums is sample-wise.
    return scale * bxdf->Sample_f(wo, wi, sample, pdf, sampledType);
}

Spectrum ScaledBxDF::rho(const Vector3f &wo, int nSamples, const Point2f *samples) const {
    // The product of 2 spectrums is sample-wise.
    return scale * bxdf->rho(wo, nSamples, samples);
}

Spectrum ScaledBxDF::rho(int nSamples, const Point2f *samples1, const Point2f *samples2) const {
    // The product of 2 spectrums is sample-wise.
    return scale * bxdf->rho(nSamples, samples1, samples2);
}

Spectrum SpecularReflection::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &sample,
    Float *pdf,
    BxDFType *sampledType
) const {
    // Compute perfect specular reflection direction about the normal. The normal vector doesn't
    // need to be known because it corresponds to the vertical axis in the reflection coordinate
    // system.
    *wi = Vector3f(-wo.x, -wo.y, wo.z);
    
    // The perfect reflection direction wi is always sampled with probability 1.
    *pdf = 1;

    // The 1/|cos(wi)| factor is meant to cancel out the |cos(wi)| factor of the integrand of
    // the scattering equation, the one that places the area differential on the plane of the
    // surface.
    return fresnel->Evaluate(CosTheta(*wi)) * R / AbsCosTheta(*wi);
}

Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW, BxDFType flags) const {
    Vector3f wi = WorldToLocal(wiW);
    Vector3f wo = WorldToLocal(woW);

    // Determine whether the incident and outgoing direction vectors are on the same or
    // opposite hemispheres.
    bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;

    Spectrum f(0.f);

    // Evaluate only the BRDFs when incident and outgoing direction vectors are on the same hemisphere.
    // Evaluate only the BTDFs when they are on opposite hemispheres.
    for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i]->MatchesFlags(flags)
            && (
                (reflect && (bxdfs[i]->type & BSDF_REFLECTION))
                || (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))
            )
        ) {
            f += bxdfs[i]->f(wo, wi);
        }
    }

    return f;
}

Spectrum BSDF::Sample_f(
    const Vector3f &woWorld,
    Vector3f *wiWorld,
    const Point2f &u,
    Float *pdf,
    BxDFType type,
    BxDFType *sampledType
) const {
    int matchingComps = NumComponents(type);
    if (matchingComps == 0) {
        *pdf = 0;
        return Spectrum(0);
    }
    // Choose one of the matching component BxDFs uniformly at random. Call it k.
    int kthMatchingComp = std::min((int) std::floor(u[0] * matchingComps), matchingComps - 1);
    BxDF * bxdf = nullptr;
    int k = kthMatchingComp;
    for (int i = 0; i < nBxDFs; ++i) {
        // Among the matching BxDFs, choose the kth one.
        if (bxdfs[i]->MatchesFlags(type) && k-- == 0) {
            bxdf = bxdfs[i];
            break;
        }
    }

    // The 2D sample u will be used for sampling the chosen BxDF, but the first of
    // its components, u[0], has already been used (for sampling the set of matching
    // BxDFs). Remap it to [0,1).
    // TODO: explain.
    Point2f uRemapped(u[0] * matchingComps - kthMatchingComp, u[1]);

    // Sample chosen BxDF. We don't care about the returned radiometric spectrum f, only
    // about wi and pdf. The BSDF's radiometric spectrum corresponding to wi is computed
    // later.
    Vector3f wi, wo = WorldToLocal(woWorld);
    *pdf = 0;
    if (sampledType) {
        *sampledType = bxdf->type;
    }
    Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
    if (*pdf == 0) {
        return 0;
    }
    *wiWorld = LocalToWorld(wi);

    // At this point, pdf stores the probability with which wi was obtained after
    // sampling the chosen BxDF. But when we chose this BxDF at random, we were actually
    // sampling the overall set of directions represented by all the matching BxDFs. So
    // wi wasn't really sampled from the chosen BxDFs's distribution, it was sampled from
    // the overall distribution of the matching BxDFs, whose PDF is the average of all the
    // PDFs.
    //
    // Compute average PDF, except in the specular case: a specular BxDF has a delta
    // distribution of directions, that is, a given wo will always be mapped to a unique wi
    // with probability 1 (pdf = 1 here already).
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1) {
        for (int i = 0; i < nBxDFs; ++i) {
            if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type)) {
                *pdf += bxdfs[i]->Pdf(wo, wi);
            }
        }
    }
    if (matchingComps > 1) {
        *pdf /= matchingComps;
    }

    // Compute value of BSDF for sampled direction.
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1) {
        // As explained earlier, we don't care about f's current value, because it corresponds
        // to only a single BxDF. We want the combined spectrum of all the matching BxDFs.
        bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
        f = 0.;
        for (int i = 0; i < nBxDFs; ++i) {
            if (bxdfs[i]->MatchesFlags(type) &&
                ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION)))) {
                f += bxdfs[i]->f(wo, wi);
            }
        }
    }

    return f;
}

Spectrum BSDF::rho(const Vector3f &wo, int nSamples, const Point2f *samples, BxDFType flags) const {
    Spectrum r(0.f);
    for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i]->MatchesFlags(flags)) {
            r += bxdfs[i]->rho(wo, nSamples, samples);
        }
    }
    return r;
}

Spectrum BSDF::rho(int nSamples, const Point2f *samples1, const Point2f *samples2, BxDFType flags ) const {
    Spectrum r(0.f);
    for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i]->MatchesFlags(flags)) {
            r += bxdfs[i]->rho(nSamples, samples1, samples2);
        }
    }
    return r;
}

Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld, BxDFType flags) const {
    if (nBxDFs == 0) {
        return 0.f;
    }

    Vector3f wo = WorldToLocal(woWorld);
    Vector3f wi = WorldToLocal(wiWorld);
    if (wo.z == 0) {
        return 0.f;
    }

    Float pdf = 0.f;
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i) {
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
    }

    // The probability of sampling wi for a given wo is the average of the PDFs of all 
    // the BxDFs that match the input flags.
    return matchingComps > 0 ? pdf / matchingComps : 0.f;
}