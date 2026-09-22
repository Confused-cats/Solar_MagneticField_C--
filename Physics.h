#pragma once
#include <vector>
#include <glm/glm.hpp>

enum class Particle_Type {
    CosmicRayParticle,
    ElectromagneticParticle
};

struct PlasmaParticle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double charge;
    double mass;
    glm::dvec3 color;
    double lifetime;
    bool active;
    Particle_Type type;
    bool escaped;

    bool gcaInitialized = false;
    double mu = 0.0;
    double v_parallel = 0.0;
};

inline constexpr double G = 6.67430e-11; //Gravittional Constant
inline constexpr double c = 299792458.0; //speed of light in m/s
inline constexpr double mu_0 = 4.0 * 3.14159265358979323846 * 1e-7; // space permeability or smthn

inline constexpr double STAR_MASS = 1.989e30; //Mass of star    
inline constexpr double STAR_RADIUS = 6.957e8; //pulsar = 10000 (10 km) sun = 696.0e9 (696000 km)
inline constexpr double SIM_LIMIT = 2.05e11; //Simulation Limit in meters
inline constexpr double STAR_OMEGA = 2.86e-6; //angular velocity (sun = 2.86e-6)
inline constexpr double globalStrength = 1e-3; // Surface Magnetic field strength 

//Atomic Masses in Kg
inline constexpr double PROTON_MASS = 1.672e-27;
inline constexpr double ELECTRON_MASS = PROTON_MASS / 1836.0; 
inline constexpr double NEUTRON_MASS = 1.6749e-27;

//... Common sense
inline constexpr double Proton_Charge = 1.602e-19;
inline constexpr double Electron_Charge = -1.602e-19;
inline constexpr double Neutron_Charge  = 0.0;

inline constexpr double alpha = 1.265; // angle of the magnetic poles (in radians)

inline constexpr double MASS_ACCRETION_RATE = 1.0e5; //experimental


/*inline double dopple_effect(double velocity) {
    velocity = 
}*/

inline constexpr double g_pull(double r) {
    if (r < STAR_RADIUS) return 0.0;
    return G * STAR_MASS / (r*r);
}

inline constexpr double r_s = (2.0 * G * STAR_MASS) / (c * c); // holder, just in case

int identifySourceSunspot(const PlasmaParticle& p);

inline double getStepSize(const glm::dvec3& pos) {
    double scale = glm::length(pos);
    return std::max(scale * 1e-5, 100.0);
}

template<typename scalarFieldFunction>
inline double scalarLaplacian(const glm::dvec3& pos, scalarFieldFunction field) {
    double h = getStepSize(pos);
    double h_sq = h * h;
    double center = field(pos);

    double d2x = (field(glm::dvec3(pos.x + h, pos.y, pos.z)) - 2.0 * center + field(glm::dvec3(pos.x - h, pos.y, pos.z))) / h_sq;
    double d2y = (field(glm::dvec3(pos.x, pos.y + h, pos.z)) - 2.0 * center + field(glm::dvec3(pos.x, pos.y - h, pos.z))) / h_sq;
    double d2z = (field(glm::dvec3(pos.x, pos.y, pos.z + h)) - 2.0 * center + field(glm::dvec3(pos.x, pos.y, pos.z - h))) / h_sq;

    return d2x + d2y + d2z;
}

template<typename vectorFieldFunction> 
inline glm::dvec3 vectorLaplacian(const glm::dvec3& pos, vectorFieldFunction field) {
    double h = getStepSize(pos);
    double h_sq = h * h;
    glm::dvec3 center = field(pos);

    glm::dvec3 d2x = (field(glm::dvec3(pos.x + h, pos.y, pos.z)) - 2.0 * center + field(glm::dvec3(pos.x - h, pos.y, pos.z))) / h_sq;
    glm::dvec3 d2y = (field(glm::dvec3(pos.x, pos.y + h, pos.z)) - 2.0 * center + field(glm::dvec3(pos.x, pos.y - h, pos.z))) / h_sq;
    glm::dvec3 d2z = (field(glm::dvec3(pos.x, pos.y, pos.z + h)) - 2.0 * center + field(glm::dvec3(pos.x, pos.y, pos.z - h))) / h_sq;

    return d2x + d2y + d2z;
}

template<typename scalarFieldFunction>
inline glm::dvec3 gradient(const glm::dvec3& pos, scalarFieldFunction field) {
    double h = getStepSize(pos);
    double two_h = 2.0 * h;

    double dx = (field(glm::dvec3(pos.x + h, pos.y, pos.z)) - field(glm::dvec3(pos.x - h, pos.y, pos.z))) / two_h;
    double dy = (field(glm::dvec3(pos.x, pos.y + h, pos.z)) - field(glm::dvec3(pos.x, pos.y - h, pos.z))) / two_h;
    double dz = (field(glm::dvec3(pos.x, pos.y, pos.z + h)) - field(glm::dvec3(pos.x, pos.y, pos.z - h))) / two_h;

    return glm::dvec3(dx, dy, dz);
}

template<typename vectorFieldFunction>  
inline double divergence(const glm::dvec3& pos, vectorFieldFunction field) {
    double h = getStepSize(pos);
    double two_h = 2.0 * h;

    double dx = (field(glm::dvec3(pos.x + h, pos.y, pos.z)).x - field(glm::dvec3(pos.x - h, pos.y, pos.z)).x) / two_h;
    double dy = (field(glm::dvec3(pos.x, pos.y + h, pos.z)).y - field(glm::dvec3(pos.x, pos.y - h, pos.z)).y) / two_h;
    double dz = (field(glm::dvec3(pos.x, pos.y, pos.z + h)).z - field(glm::dvec3(pos.x, pos.y, pos.z - h)).z) / two_h;

    return dx + dy + dz;
}

template<typename vectorFieldFunction>
inline glm::dvec3 curl(const glm::dvec3& pos, vectorFieldFunction field) {
    double h = getStepSize(pos);
    double two_h = 2.0 * h;

    double dFz_dy = (field(glm::dvec3(pos.x, pos.y + h, pos.z)).z - field(glm::dvec3(pos.x, pos.y - h, pos.z)).z) / two_h;
    double dFy_dz = (field(glm::dvec3(pos.x, pos.y, pos.z + h)).y - field(glm::dvec3(pos.x, pos.y, pos.z - h)).y) / two_h;

    double dFx_dz = (field(glm::dvec3(pos.x, pos.y, pos.z + h)).x - field(glm::dvec3(pos.x, pos.y, pos.z - h)).x) / two_h;
    double dFz_dx = (field(glm::dvec3(pos.x + h, pos.y, pos.z)).z - field(glm::dvec3(pos.x - h, pos.y, pos.z)).z) / two_h;

    double dFx_dy = (field(glm::dvec3(pos.x, pos.y + h, pos.z)).x - field(glm::dvec3(pos.x, pos.y - h, pos.z)).x) / two_h;
    double dFy_dx = (field(glm::dvec3(pos.x + h, pos.y, pos.z)).y - field(glm::dvec3(pos.x - h, pos.y, pos.z)).y) / two_h;

    return glm::dvec3(
        dFz_dy - dFy_dz,
        dFx_dz - dFz_dx,
        dFy_dx - dFx_dy
    );
}


// curl (glm::dvec3(pos.x, pos.y, pos.z)) jfewuobiusbdgdipqb
// Extern mutable state
extern double starRotationAngle;
extern bool useParkerSpiral;
extern bool enableElectricField; 
extern std::vector<PlasmaParticle> plasmaParticles; // Shared state

void updateCachedSunspotPos();

double calculateAlfvenRadius();

// Field Calculations
glm::dvec3 calculateMagneticField(const glm::dvec3& position);
glm::dvec3 calculateElectricField(const glm::dvec3& position, const glm::dvec3& B);

inline double getFieldMagnitude(const glm::dvec3& pos) {
    return glm::length(calculateMagneticField(pos));
}

glm::dvec3 calculateCurrentDensity(const glm::dvec3& position);

double calculateMagneticMoment(const glm::dvec3& velocity, double mass, const glm::dvec3& B);
glm::dvec3 calculateGradBDrift(const glm::dvec3& pos, double charge, double mu);
double calculateMirrorForce(const glm::dvec3& pos, double mu);

// Particle Emitters & Injectors
void injectPlasmaField(std::vector<PlasmaParticle>& particles);
void injectCosmicRay(std::vector<PlasmaParticle>& particles);
void triggerCoronalMassEjection(std::vector<PlasmaParticle>& particles);
void triggerCarringtonEvent(std::vector<PlasmaParticle>& particles);
void triggerLocalizedCME(std::vector<PlasmaParticle>& particles, const glm::dvec3& spotPos, double calculatedSpeed, int particleCount);
void evaluateMHDInstabilities(std::vector<PlasmaParticle>& particles, double physicsDt);

// Collision & Environment Interactions
void handleCollisions(std::vector<PlasmaParticle>& particles);
// Simulation Solvers (Boris Integrator)
void updateChargedParticle(PlasmaParticle& p, double dt, double starRadius);
void updateParticles(std::vector<PlasmaParticle>& particles, double dt);