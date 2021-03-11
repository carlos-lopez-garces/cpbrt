#include "cpbrt.h"
#include "reflection.h"
#include "sampling.h"

Spectrum BxDF::Sample_f(
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
    BxDFType *sampledType = nullptr
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