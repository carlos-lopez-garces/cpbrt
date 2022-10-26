#include "cpbrt.h"
#include "medium.h"
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
    return FrDielectric(cosThetaI, etaI, etaT);
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
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
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
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
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

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
Spectrum LayeredBxDF<TopBxDF, BottomBxDF, twoSided>::f(const Vector3f &wo, const Vector3f &wi) const {
    SampledSpectrum f(0.f);

    if(twoSided && wo.z < 0) {
        // If the surface is twoSided, the TopBxDF is always the top, regardless of the side of
        // incidence. In the shading or reflection coordinate system, both wo and wi point outward and away
        // from the surface. When wo.z < 0, we reverse the vectors so that the side of incidence of the
        // surface becomes the top layer.
        wo = -wo;
        wi = -wi;
    }

    // Determine entry and exit layers.

    TopOrBottomBxDF<TopBxDF, BottomBxDF> enterInterface;
    bool enteredTop = twoSided || wo.z > 0;
    if (enteredTop) {
        enterInterface = &top;
    } else {
        enterInterface = &bottom;
    }

    TopOrBottomBxDF<TopBxDF, BottomBxDF> exitInterface, nonExitInterface;
    // If the incident and outgoing vectors lie on the same hemisphere, the ray is being reflected.
    // If the ray is being reflected after arriving first at the bottom layer (e.g. a transmitted ray that
    // gets reflected at the exit interface, like in the total internal reflection case), then the 
    // bottom layer is both the entry and exit interfaces (and the top layer is none of them).
    // [First conditional branch executes.]
    //
    // If the ray is being reflected after arriving first at the top layer, then the top layer is both
    // the entry and exit interfaces (and the bottom layer is none of them).
    // [Second conditional branch executes.]
    //
    // If the incident and outgoing vectors don't lie on the same hemisphere, the ray is being transmitted.
    // If the ray is being transmitted after arriving first at the top layer, then the top layer is the
    // entry interface and the bottom layer is the exit interface.
    // [First conditional branch executes.]
    //
    // If the ray is being transmitted after arriving first at the bottom layer, then ... what?
    if (SameHemisphere(wo, wi) ^ enteredTop) {
        exitInterface = &bottom;
        nonExitInterface = &top;
    } else {
        exitInterface = &top;
        nonExitInterface = &bottom;
    }

    // The thickness of the interface of a layered material is not negligible and we must account
    // for the transmission of the ray through this "coat".
    //
    // If the ray is either being reflected or it's being transmitted after entering through the
    // top layer [not sure about the latter], but not both, the traveled distance through the coat
    // is actually 0; otherwise, the traveled distance is measured by the thickness of the layer.
    Float exitZ = (SameHemisphere(wo, wi) ^ enteredTop) ? 0 : thickness;

    if (SameHemisphere(wo, wi)) {
        // The ray is being reflected. Sample the top layer's BRDF.
        f = nSamples * enterInterface.f(wo, wi);
    }

    // TODO: might need to enhance the initialization of RNG.
    RNG rng;

    for (int s = 0; s < nSamples; ++s) {
        // Sample the entry layer's BTDF. Note that we sample it using wo.
        Vector3f sampledWi;
        Float woPdf = 0.f;
        Spectrum fWo = enterInterface.Sample_f(
            wo, &sampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &woPdf, BxDFType::BSDF_TRANSMISSION
        );
        if (fWo.IsBlack() || woPdf == 0.f || sampledWi.z == 0) {
            // TODO: when would the sampled incident directon wis's z coordinate be 0?
            continue;
        }

        // Sample the exit layer's BTDF. Note that we sample it using wi.
        Vector3f sampledWo;
        Float wiPdf = 0.f;
        Spectrum fWi = exitInterface.Sample_f(
            wi, &sampledWo, Point2f(rng.UniformFloat(), rng.UniformFloat()), &wiPdf, BxDFType::BSDF_TRANSMISSION
        );
        if (fWi.IsBlack() || wiPdf == 0.f || sampledWo.z == 0) {
            // TODO: when would the z coordinate be 0?
            continue;
        }

        Spectrum beta = fWo * AbsCosTheta(sampledWi) / woPdf;
        Float z = enteredTop ? thickness : 0;
        Vector3f w = sampledWi;
        HenyeyGreensteinPhaseFunction phase(g);

        // Random walk.
        for (int depth = 0; depth < maxDepth; ++depth) {
            if (depth > 3 && beta.MaxComponentValue() < 0.25f) {
                Float q = std::max<Float>(0, 1 - beta.MaxComponentValue());

                if (rng.UniformFloat() < q) {
                    // Terminate random walk using Russian roulette.
                    break;
                }

                beta /= 1 - q;
            }

            // Transmission through the medium between the layers.
            if (!albedo) {
                // No absorption means that the radiance carried by the ray through the exit layer's
                // interface will be the same it carried when it entered the entry layer's interface.
                //
                // Advance z to the next layer's interface.
                z = (z == thickness) ? 0 : thickness;
                beta *= Tr(thickness, w);
            } else {
                // Sample a potential scattering event.

                Float sigma_t = 1;

                // Sample a exponential random variable with parameter lambda = sigma_t / |w.z|. 
                Float dz = SampleExponential(rng.UniformFloat(), sigma_t / std::abs(w.z));

                Float zp = w.z > 0 ? (z + dz) : (z - dz);
                if (zp == z) {
                    return Spectrum(0);
                }

                if (0 < zp && zp < thickness) {
                    // Process scattering event.
                    Float wt = 1;
                    if (!(exitInterface.type & BxDFType::BSDF_SPECULAR)) {
                        wt = PowerHeuristic(1, wiPdf, 1, phase.Pdf(-w, -sampledWo));
                    }

                    f += beta * albedo * phase.p(-w, -sampledWo) * wt * Tr(zp - exitZ, sampledWo) * fWi / wiPdf;

                    // Sample phase function.
                    Vector3f phaseWi;
                    Float phaseSample = phase.Sample_p(-w, &phaseWi, Point2f(rng.UniformFloat(), rng.UniformFloat()));
                    if (phaseWi.z == 0) {
                        continue;
                    }
                    beta *= albedo * phaseSample / phase.Pdf(-w, phaseWi);
                    w = phaseWi;
                    z = zp;

                    // Process scattering through exit interface.
                    if (((z < exitZ && w.z > 0) || (z > exitZ && w.z < 0)) && !(exitInterface.type & BxDFType::BSDF_SPECULAR)) {
                        Spectrum fExit = exitInterface.f(-w, wi);
                        if (fExit) {
                            Float exitPDF = exitInterface.Pdf(-w, wi);
                            Float wt = PowerHeuristic(1, phase.p(-w, -sampledWo), 1, exitPDF);
                            f += beta * Tr(zp - exitZ, phaseWi) * fExit * wt;
                        }
                    }

                    continue;
                }

                z = Clamp(zp, 0, thickness);
            }

            if (z == exitZ) {
                // Process scattering at exit interface.
                Vector3f sampledExitWi;
                Float exitPdf = 0.f;
                Spectrum fExit = exitInterface.Sample_f(
                    -w, &sampledExitWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &exitPdf, BxDFType::BSDF_REFLECTION
                );
                if (fExit.IsBlack() || exitPdf == 0.f || sampledExitWi.z == 0) {
                    break;
                }

                beta *= fExit * AbsCosTheta(sampledExitWi) / exitPdf;
                w = sampledExitWi;
            } else {
                // Process scattering at nonexit interface.
                if (!(nonExitInterface.type & BxDFType::BSDF_SPECULAR)) {
                    Float wt = 1;
                    if (!(exitInterface.type & BxDFType::BSDF_SPECULAR)) {
                        wt = PowerHeuristic(1, wiPdf, 1, nonExitInterface.Pdf(-w, -sampledWo));
                    }

                    f += beta * nonExitInterface.f(-w, -sampledWo) * AbsCosTheta(sampledWo) * wt * Tr(thickness, sampledWo) * fWi / wiPdf;
                }

                // Sample new direction using BSDF at nonexit interface.
                Vector3f sampledNonExitWo;
                Float nonExitPdf = 0.f;
                Spectrum fNonExit = exitInterface.Sample_f(
                    -w, &sampledNonExitWo, Point2f(rng.UniformFloat(), rng.UniformFloat()), &nonExitPdf, BxDFType::BSDF_REFLECTION
                );
                if (fNonExit.IsBlack() || nonExitPdf == 0.f || sampledNonExitWo.z == 0) {
                    break;
                }

                beta *= fNonExit * AbsCosTheta(sampledNonExitWo) / nonExitPdf;
                w = sampledNonExitWo;

                if (!(exitInterface.type & BxDFType::BSDF_SPECULAR)) {
                    Spectrum fExit = exitInterface.f(-w, wi);
                    if (!fExit.IsBlack()) {
                        Float wt = 1;
                        if (!(nonExitInterface.type & BxDFType::BSDF_SPECULAR)) {
                            Float exitPdf = exitInterface.Pdf(-w, wi);
                            wt = PowerHeuristic(1, nonExitPdf, 1, exitPdf);
                        }

                        f += beta * Tr(thickness, sampledNonExitWo) * fExit * wt;
                    }
                }
            }
        }
    }

    return f / nSamples;
}

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
Spectrum LayeredBxDF<TopBxDF, BottomBxDF, twoSided>::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &sample,
    Float *pdf,
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
) const {
    bool flipWi = false;
    if(twoSided && wo.z < 0) {
        // If the surface is twoSided, the TopBxDF is always the top, regardless of the side of
        // incidence. In the shading or reflection coordinate system, both wo and wi point outward and away
        // from the surface. When wo.z < 0, we reverse the vectors so that the side of incidence of the
        // surface becomes the top layer.
        wo = -wo;
        flipWi = true;
    }

    bool enteredTop = twoSided || wo.z > 0;
    Vector3f sampledWi;
    Float entryPdf = 0.f;
    BxDFType entrySampledType;
    Spectrum entryF = enteredTop 
        ? top.Sample_f(wo, &sampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &entryPdf, &entrySampledType, uc, mode)
        : bottom.Sample_f(wo, &sampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &entryPdf, &entrySampledType, uc, mode);
    if (entryF.IsBlack() || entryPdf == 0.f || sampledWi.z == 0) {
        return Spectrum(0.f);
    }
    if (entrySampledType & BSDF_REFLECTION) {
        if (flipWi) {
            sampledWi = -sampledWi;
        }
        *wi = sampledWi;
        if (sampledType) {
            *sampledType = entrySampledType;
        }
        return entryF;
    }

    // TODO: might need to enhance the initialization of RNG.
    RNG rng;

    Vector3f retW = sampledWi;
    bool specularPath = entrySampledType & BSDF_SPECULAR;
    Spectrum retF = entryF * AbsCosTheta(sampledWi);
    Float retPdf = entryPdf;
    Float z = enteredTop ? thickness : 0;
    HenyeyGreensteinPhaseFunction phase(g);

    // Random walk.
    for (int depth = 0; depth < maxDepth; ++depth) {
        Float rrBeta = retF.MaxComponentValue() / retPdf;
        if (depth > 3 && rrBeta < 0.25f) {
            Float q = std::max<Float>(0, 1 - rrBeta);

            if (rng.UniformFloat() < q) {
                // Terminate random walk using Russian roulette.
                return Spectrum(0.f);
            }

            retPdf *= 1 - q;
        }

        if (retW.z == 0) {
            return Spectrum(0.f);
        }

        // Transmission through the medium between the layers.
        if (!albedo) {
            // No absorption means that the radiance carried by the ray through the exit layer's
            // interface will be the same it carried when it entered the entry layer's interface.
            //
            // Advance z to the next layer's interface.
            z = (z == thickness) ? 0 : thickness;
            retF *= Tr(thickness, retW);
        } else {
            // Sample a potential scattering event.

            Float sigma_t = 1;

            // Sample a exponential random variable with parameter lambda = sigma_t / |w.z|. 
            Float dz = SampleExponential(rng.UniformFloat(), sigma_t / AbsCosTheta(retW));

            Float zp = retW.z > 0 ? (z + dz) : (z - dz);
            if (zp == z) {
                return Spectrum(0);
            }

            if (0 < zp && zp < thickness) {
                // Process scattering event.
                
                // Sample phase function.
                Vector3f phaseWi;
                Float phaseSample = phase.Sample_p(-retW, &phaseWi, Point2f(rng.UniformFloat(), rng.UniformFloat()));
                Float phasePdf = phase.Pdf(-retW, phaseWi);
                if (phasePdf == 0.f || phaseWi.z == 0) {
                    return Spectrum(0.f);
                }
                retF *= albedo * phaseSample;
                retPdf *= phasePdf;
                specularPath = false;
                retW = phaseWi;
                z = zp;

                continue;
            }

            z = Clamp(zp, 0, thickness);
        }

        TopOrBottomBxDF<TopBxDF, BottomBxDF> interface;
        if (z == 0) {
            interface = &bottom;
        } else {
            interface = &top;
        }

        Vector3f interfaceSampledWi;
        Float interfacePdf = 0.f;
        BxDFType interfaceSampledType;
        Spectrum interfaceF = interface.Sample_f(-retW, &interfaceSampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &interfacePdf, &interfaceSampledType, uc, mode);
        if (interfaceF.IsBlack() || interfacePdf == 0.f || interfaceSampledWi.z == 0) {
            return Spectrum(0.f);
        }
        retF *= interfaceF;
        retPdf *= interfacePdf;
        specularPath &= interfaceSampledType & BSDF_SPECULAR;
        retW = interfaceSampledWi;

        if (interfaceSampledType & BSDF_TRANSMISSION) {
            BxDFType type = SameHemisphere(wo, retW) ? BSDF_REFLECTION : BSDF_TRANSMISSION;
            type |= specularPath ? BSDF_SPECULAR : BSDF_GLOSSY;
            if (flipWi) {
                retW = -retW;
            }
            *wi = retW;
            *pdf = retPdf;
            *sampledType = interfaceSampledType; 
            return retF;
        }

        retF *= AbsCosTheta(interfaceSampledWi);
    }

    return Spectrum(0.f);
}

template <typename TopBxDF, typename BottomBxDF, bool twoSided>
Float LayeredBxDF<TopBxDF, BottomBxDF, twoSided>::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (twoSided && wo.z < 0) {
        wo = -wo;
        wi = -wi;
    }

    RNG rng;

    bool enteredTop = twoSided || wo.z > 0;
    Float pdfSum = 0.f;
    if (SameHemisphere(wo, wi)) {
        pdfSum += enteredTop ? nSamples * top.Pdf(wo, wi) : nSamples * bottom.Pdf(wo, wi);
    }

    for (int s = 0; s < nSamples; ++s) {
        if (SameHemisphere(wo, wi)) {
            TopOrBottomBxDF<TopBxDF, BottomBxDF> rInterface;
            TopOrBottomBxDF<TopBxDF, BottomBxDF> tInterface;

            if (enteredTop) {
                rInterface = &bottom;
                tInterface = &top;
            } else {
                rInterface = &top;
                tInterface = &bottom;
            }

            BxDFType trans = BSDF_TRANSMISSION;

            Vector3f sampledWi;
            Float pdfWo;
            Spectrum fWo = tInterface.Sample_f(wo, &sampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &pdfWo);

            Vector3f sampledWo;
            Float pdfWi;
            Spectrum fWi = tInterface.Sample_f(wi, &sampledWo, Point2f(rng.UniformFloat(), rng.UniformFloat()), &pdfWi);

            if (!fWo.IsBlack() && pdfWo > 0 && !fWi.IsBlack() && pdfWi > 0) {
                // If nonspecular.
                if (tInterface.type & (BSDF_DIFFUSE | BSDF_GLOSSY)) {
                    pdfSum += rInterface.Pdf(-sampledWi, -sampledWo);
                } else {
                    // MIS for product.
                    Vector3f rsSampledWi;
                    Float pdfRS;
                    Spectrum fRS = rInterface.Sample_f(-sampledWi, &rsSampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &pdfRS);
                    if (!fRS.IsBlack() && pdfRS > 0) {
                        if (!(rInterface.type & (BSDF_DIFFUSE | BSDF_GLOSSY))) {
                            pdfSum += tInterface.Pdf(-rsSampledWi, wi);
                        } else {
                            Float rPdf = rInterface.Pdf(-sampledWi, -sampledWo);
                            Float wt = PowerHeuristic(1, pdfWi, 1, rPdf);
                            pdfSum += wt * rPdf;

                            Float tPdf = tInterface.Pdf(-rsSampledWi, wi);
                            wt = PowerHeuristic(1, pdfRS, 1, tPdf);
                            pdfSum += wt * tPdf;
                        }
                    }
                }
            }
        } else {
            TopOrBottomBxDF<TopBxDF, BottomBxDF> toInterface;
            TopOrBottomBxDF<TopBxDF, BottomBxDF> tiInterface;
            if (enteredTop) {
                toInterface = &top;
                tiInterface = &bottom;
            } else {
                toInterface = &bottom;
                tiInterface = &top;
            }

            Vector3f toSampledWi;
            Float toPdf;
            BxDFType toSampledType;
            Spectrum fTo = toInterface.Sample_f(wo, &toSampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &toPdf, &toSampledType, rng.UniformFloat());
            if (fTo.IsBlack() || toPdf == 0 || toSampledWi.z == 0 || toSampledType & BSDF_REFLECTION) {
                continue;
            }

            Vector3f tiSampledWi;
            Float tiPdf;
            BxDFType tiSampledType;
            Spectrum fTi = tiInterface.Sample_f(wi, &tiSampledWi, Point2f(rng.UniformFloat(), rng.UniformFloat()), &tiPdf, &tiSampledType, rng.UniformFloat());
            if (fTi.IsBlack() || tiPdf == 0 || tiSampledWi.z == 0 || tiSampledType & BSDF_REFLECTION) {
                continue;
            }

            if (toInterface.type & BSDF_SPECULAR) {
                pdfSum += tiInterface.Pdf(-toSampledWi, wi);
            } else if (tiInterface.type & BSDF_SPECULAR) {
                pdfSum += toInterface.Pdf(wo, -tiSampledWi);
            } else {
                pdfSum += (toInterface.Pdf(wo, -tiSampledWi) + tiInterface.Pdf(-toSampledWi, wi)) / 2;
            }
        }
    }

    return Lerp(0.9f, 1 / (4 * Pi), pdfSum / nSamples);
}

Spectrum SpecularReflection::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &sample,
    Float *pdf,
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
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

Spectrum SpecularTransmission::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &sample,
    Float *pdf,
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
) const {
    // Is the ray entering or exiting the medium? cos(theta) is computed in reflection space;
    // theta is measured from the transmission boundary's normal; in reflection space, cos(theta)
    // corresponds to the z coordinate of wo.
    bool entering = CosTheta(wo) > 0;
    Float etaI = entering ? etaA : etaB;
    Float etaT = entering ? etaB : etaA;

    // Compute ray direction for specular transmission.
    if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi)) {
        return 0;
    }

    // The distribution of transmission directions is a Dirac delta. For a given transmission
    // direction, a unique incident direction is determined by Snell's law deterministically.
    // This direction is "sampled" with probability 1.
    *pdf = 1;

    // T * (1 - Fr), where T is a scaling factor. The 1 - Fr is larger for incidence directions
    // near the normal. Transmisison is thus stronger when the ray enters the medium orthogonally
    // than when it does at a grazing angle. The stronger transmission, the more visible the
    // transmission medium becomes (or what's on the other side of the boundary); the weaker
    // transmission is, the stronger reflection becomes.
    Spectrum ft = T * (Spectrum(1.) - fresnel.Evaluate(CosTheta(*wi)));

    // Account for non-symmetry with transmission to different medium.
    // TODO: explain after reading section 16.1.3.
    if (mode == TransportMode::Radiance) {
        ft *= (etaI * etaI) / (etaT * etaT);
    }

    // The 1/|cos(wi)| factor is meant to cancel out the |cos(wi)| factor of the integrand of
    // the scattering equation, the one that places the area differential on the plane of the
    // transmission medium's boundary.
    return ft / AbsCosTheta(*wi);
}

Spectrum FresnelSpecularReflectionTransmission::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &u,
    Float *pdf,
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
) const {
    // Evaluate the Fresnel reflectance equation between 2 dielectric media. This Fresnel term
    // is used to "modulate" the contributions of the BRDF and the BTDF; in reality, F is
    // interpreted as a probability of sampling the BRDF vs the BTDF in a given call. For example,
    // at glancing angles, where reflection is strong, the BRDF has a larger contribution, so it
    // has a larger probability of being sampled.
    Float F = FrDielectric(CosTheta(wo), etaA, etaB);
    // u is a sample from a uniform distribution that is used as a probability threshold to choose
    // between sampling the BRDF (< u) and the BTDF(> u).
    if (u[0] < F) {
        // Evaluate the BRDF.

        // Perfect specular reflection direction.
        *wi = Vector3f(-wo.x, -wo.y, wo.z);
        if (sampledType) {
            *sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
        }

        *pdf = F;

        // Same as SpecularReflection::Sample_f.
        return F * R / AbsCosTheta(*wi);
    } else {
        // Evaluate the BTDF. 
        // Is the ray entering or exiting the medium? cos(theta) is computed in reflection space;
        // theta is measured from the transmission boundary's normal; in reflection space, cos(theta)
        // corresponds to the z coordinate of wo.
        bool entering = CosTheta(wo) > 0;
        Float etaI = entering ? etaA : etaB;
        Float etaT = entering ? etaB : etaA;

        // Perfect specular transmission direction.
        if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi)) {
            return 0;
        }

        // T * (1 - Fr), where T is a scaling factor. The 1 - Fr is larger for incidence directions
        // near the normal. Transmisison is thus stronger when the ray enters the medium orthogonally
        // than when it does at a grazing angle. The stronger transmission, the more visible the
        // transmission medium becomes (or what's on the other side of the boundary); the weaker
        // transmission is, the stronger reflection becomes.
        Spectrum ft = T * (1 - F);

        if (sampledType) {
            *sampledType = BxDFType(BSDF_SPECULAR | BSDF_TRANSMISSION);
        }

        *pdf = 1 - F;

        // Same as SpecularTransmission::Sample_f.
        return ft / AbsCosTheta(*wi);
    }
}

Spectrum OrenNayarReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    // ThetaI (ThetaO) is the colatitude angle of the spherical coordinate of wi (wo) in the
    // shading coordinate system. See reflection.h.
    Float sinThetaI = SinTheta(wi);
    Float sinThetaO = SinTheta(wo);

    // Compute cosine term of Oren-Nayar model.
    Float maxCos = 0;
    if (sinThetaI > 1e-4 && sinThetaO > 1e-4) {
        // PhiI (PhiO) is the polar angle of the spherical coordinate of wi (wo).
        Float sinPhiI = SinPhi(wi);
        Float cosPhiI = CosPhi(wi);
        Float sinPhiO = SinPhi(wo);
        Float cosPhiO = CosPhi(wo);

        // cos(PhiI - PhiO) computed using the trigonometric identity cos(PhiI)cos(PhiO) + sin(PhiI)sin(PhiO).
        Float dCos = cosPhiI*cosPhiO + sinPhiI*sinPhiO;
        maxCos = std::max((Float) 0, dCos);
    }

    // Compute sine and tangent terms of Oren-Nayar model.
    // alpha = max(ThetaI, ThetaO) = min(cos(ThetaI), cos(ThetaO))
    // beta = min(ThetaI, ThetaO) = max(cos(ThetaI), cos(ThetaO)).
    Float sinAlpha;
    Float tanBeta;
    if (AbsCosTheta(wi) > AbsCosTheta(wo)) {
        sinAlpha = sinThetaO;
        tanBeta = sinThetaI / AbsCosTheta(wi);
    } else {
        sinAlpha = sinThetaI;
        tanBeta = sinThetaO / AbsCosTheta(wo);
    }

    return R * InvPi * (A + B*maxCos*sinAlpha*tanBeta);
}

Spectrum TorranceSparrowMicrofacetReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);

    // Half-angle.
    Vector3f wh = wi + wo;

    // Edge cases: perfectly grazing incident or outgoing directions.
    if (cosThetaI == 0 || cosThetaO == 0) {
        return Spectrum(0.);
    }
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) {
        return Spectrum(0.);
    }

    wh = Normalize(wh);

    // Fresnel reflectance; fraction of incident light that gets reflected in the half-angle direction.
    Spectrum F = fresnel->Evaluate(Dot(wi, Faceforward(wh, Vector3f(0,0,1))));

    // Torrance-Sparrow BRDF, where D(wh) gives the fraction of differential area covered by microfacets
    // with normal equal to the half-angle (only those reflect light for given wo and wi) and G(wo, wi)
    // is the geometric attenuation factor that accounts for masking and shadowing. 
    return R * distribution->D(wh) * distribution->G(wo, wi) * F / (4 * cosThetaI * cosThetaO);
}

// TODO.
Spectrum TorranceSparrowMicrofacetReflection::Sample_f(
    const Vector3f &wo,
    Vector3f *wi,
    const Point2f &u,
    Float *pdf,
    BxDFType *sampledType,
    Float uc,
    TransportMode mode
) const {
    // Sample microfacet orientation $\wh$ and reflected direction $\wi$
    if (wo.z == 0) return 0.;
    Vector3f wh = distribution->Sample_wh(wo, u);
    if (Dot(wo, wh) < 0) return 0.;   // Should be rare
    *wi = Reflect(wo, wh);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

    // Compute PDF of _wi_ for microfacet reflection
    *pdf = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
    return f(wo, *wi);
}

// TODO.
Float TorranceSparrowMicrofacetReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;
    Vector3f wh = Normalize(wo + wi);
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
}

// The Torrance-Sparrow BTDF.
Spectrum TorranceSparrowMicrofacetTransmission::f(const Vector3f &wo, const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) {
        return Spectrum(0);
    }

    Float cosThetaO = CosTheta(wo);
    Float cosThetaI = CosTheta(wi);
    // Theta in the spherical coordinate (Theta, Phi, r) is measured from the normal of the surface
    // at the shaded point (in shading space) to vector wo or vector wi. These vectors are normalized.
    // cos(Theta) is the scalar component of wo or wi in the direction of the normal and corresponds
    // to wo.z or wi.z because they are normalized.
    if (cosThetaO == 0 || cosThetaI == 0) {
        // wo or wi is tangential to the surface. No transmission.
        return Spectrum(0);
    }

    Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);

    // Half vector.
    Vector3f wh = Normalize(wo + wi * eta);
    if (wh.z < 0) {
        wh = -wh;
    }

    if (Dot(wo, wh)*Dot(wi, wh) > 0) {
        // wo, wi, and wh are all on the same side hemisphere. This can't be because wo and wi
        // are on different hemispheres during transmission.
        return Spectrum(0);
    }

    Spectrum F = fresnel.Evaluate(Dot(wo, wh));
    Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
    Float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;

    // BTDF.
    return (Spectrum(1.f) - F) * T * std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta * AbsDot(wi, wh) * AbsDot(wo, wh) * factor * factor / (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}

Spectrum TorranceSparrowMicrofacetTransmission::Sample_f(
    const Vector3f &wo, Vector3f *wi, const Point2f &u, Float *pdf, BxDFType *sampledType
) const {
    if (wo.z == 0) {
        return 0.;
    }

    Vector3f wh = distribution->Sample_wh(wo, u);
    if (Dot(wo, wh) < 0) {
        return 0.;
    }

    Float eta = CosTheta(wo) > 0 ? (etaA / etaB) : (etaB / etaA);
    if (!Refract(wo, (Normal3f)wh, eta, wi)) {
        return 0;
    }

    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float TorranceSparrowMicrofacetTransmission::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) {
        return 0;
    }
    Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = Normalize(wo + wi * eta);

    if (Dot(wo, wh) * Dot(wi, wh) > 0) {
        return 0;
    }

    Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
    Float dwh_dwi = std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));
    return distribution->Pdf(wo, wh) * dwh_dwi;
}

WardReflection::WardReflection(const Spectrum &Rd, const Spectrum &Rs, const Float Ax, const Float Ay)
	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), Rd(Rd), Rs(Rs) {
	ax = Ax;
	ay = Ay;
	invax2 = 1.f/(ax*ax);
	invay2 = 1.f/(ay*ay);
	const2 = (4*Pi*ax*ay);
}

Spectrum WardReflection::f(const Vector3f &wo, const Vector3f &wi) const{
	Spectrum specular(0.f);
	Vector3f H = wi + wo;
	if(H.z <= 0.f) return specular;
	Float const1 = wi.z*wo.z;
	if(const1 <= 0.f) return specular;
	const1 = sqrtf(const1);
	Float const3 = exp(-1.f * (H.x*H.x*invax2 + H.y*H.y*invay2)/(H.z*H.z) );
	specular = Rs * const3/(const1*const2); 
	return specular;
}

Spectrum WardReflection::f(const Vector3f &wo, const Vector3f &wi, const Vector3f &H) const{
	Spectrum specular(0.f);
	
	if(H.z <= 0.f) return specular;
	Float const1 = wi.z*wo.z;
	if(const1 <= 0.f) return specular;
	const1 = sqrtf(const1);
	Float const3 = exp(-1.f * (H.x*H.x*invax2 + H.y*H.y*invay2)/(H.z*H.z) );
	specular = Rs * const3/(const1*const2); 
	return specular;	
}

Spectrum WardReflection::Sample_f(
    const Vector3f &wo, Vector3f *wi, const Point2f &u, Float *pdf, BxDFType *sampledType
) const {
	Vector3f h;
	Float phi = atanf(ay*tan(2*Pi*u[2])/ax);
	Float cosPhi = cosf(phi);
	Float sinPhi = sqrtf(1-cosPhi*cosPhi);
	Float theta = atanf(sqrtf(-logf(u[1])/(cosPhi*cosPhi*invax2 + sinPhi*sinPhi*invay2)));
	h.z = cosf(theta);
	Float cosTheta2 = h.z*h.z;
	Float sinTheta = sqrtf(1-cosTheta2);
	Float tanTheta2 = (1-cosTheta2)/cosTheta2;
	h.x = cosPhi*sinTheta;
	h.y = sinPhi*sinTheta;
	if(Dot(wo, h) < 0.f) h = -h;
	*wi = -wo + 2.f * Dot(wo, h) * h;
	*pdf = expf(-tanTheta2*(cosPhi*cosPhi*invax2 + sinPhi*sinPhi*invay2))/(const2*Dot(h,wo)*cosTheta2*h.z);
	return f(wo, *wi, h);
}

Float WardReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
	Vector3f h = Normalize(wi+wo);
	Float cosTheta2 = h.z*h.z;
	Float tanTheta2 = (1-cosTheta2)/cosTheta2;
	Float cosPhi = CosPhi(h);
	Float sinPhi = SinPhi(h);

	return expf(-tanTheta2*(cosPhi*cosPhi*invax2 + sinPhi*sinPhi*invay2))/(const2*Dot(h,wo)*cosTheta2*h.z);
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