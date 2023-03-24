#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        reflect(wo, wi);
        *pdf = 1.0;
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        double theta_h = getTheta(h);
        return  (exp(-pow(tan(theta_h),2)/pow(alpha, 2)))/ (PI * pow(alpha, 2) * pow(cos(theta_h), 4));
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        
        return Vector3D();
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        Vector3D h = (wo + wi).unit();
        Vector3D n = Vector3D(0,0,1);

//        if (dot(n, wo) < 0 || dot(n, wi) < 0){
//            return Vector3D();
//        }

        return (F(wi) * G(wo, wi) * D(h)) / (4.0 * dot(n, wo) * dot(n, wi));

    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        *wi = cosineHemisphereSampler.get_sample(pdf);
        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        double n = wo.z < 0 ? ior : 1.0 / ior;
        bool result = refract(wo, wi, ior);
        *pdf = 1.0;
        if (!result) {
            return Vector3D();
        }
        return transmittance / abs_cos_theta(*wi) / (n * n);
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        if (!refract(wo, wi, ior)) {
            // If there is total internal reflection, wo does not have a valid refracted wi
            reflect(wo, wi);
            *pdf = 1.0;
            return reflectance / abs_cos_theta(*wi);
        }
        // Else
        double r0 = (ior - 1.0) * (ior - 1.0) / (ior + 1.0) / (ior + 1.0);
        double term = (1.0 - abs_cos_theta(*wi)) * (1.0 - abs_cos_theta(*wi)) * (1.0 - abs_cos_theta(*wi)) * (1.0 - abs_cos_theta(*wi)) * (1.0 - abs_cos_theta(*wi));
        double R = r0 + (1.0 - r0) * term;
        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            refract(wo, wi, ior);
            *pdf = 1.0 - R;
            double n = wo.z < 0 ? ior : 1.0 / ior;
            return (1.0 - R) * transmittance / abs_cos_theta(*wi) / (n*n);
        }
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.

        (*wi) = Vector3D(-wo.x, -wo.y, wo.z);
    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.

        // Check if entering or exiting
        double n, sign;
        if (wo.z < 0) {
            // Inside, exiting
            n = ior;
            sign = 1.0;
        }
        else {
            // Outside, entering
            n = 1.0 / ior;
            sign = -1.0;
        }
        double term = 1.0 - (n * n * (1.0 - (wo.z * wo.z)));
        if (term < 0) {
            return false;
        }
        *wi = Vector3D(-n * wo.x, -n * wo.y, sign * sqrt(term));
        return true;

    }

} // namespace CGL
