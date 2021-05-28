#include "cpbrt.h"
#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "scene.h"
#include "interaction.h"

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
        return f(woWorld, wiWorld, type);
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