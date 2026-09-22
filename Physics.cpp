#include "Physics.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <random>
#include <limits>

double starRotationAngle = 0.0;
bool useParkerSpiral = true;
bool enableElectricField = true;

// Cache variables, preventing redundant trig calcs
double cachedCosRotation = 1.0;
double cachedSinRotation = 0.0;
double cachedCosInverseRotation = 1.0;
double cachedSinInverseRotation = 0.0;

std::vector<PlasmaParticle> plasmaParticles;
glm::dvec3 cachedSunspotPos = glm::dvec3(STAR_RADIUS, 0.0, 0.0);

struct SunspotAnchor {
    glm::dvec3 staticPos;         
    glm::dmat3 toStaticTransform;  // Matrix mapping local (X-up) frame to star static frame
    glm::dmat3 toLocalTransform;   // Inverse matrix mapping star static frame back to local
    double tiltAngle = 0.0;

    double freeMagneticEnergy = 0.0;
    double trappedPlasmaMass = 0.0;

    double inflation = 1.0;
    double eruptionFactor =0.0;
};

std::vector<SunspotAnchor> activeSunspots;

void initializeRandomSunspots(int count) {
    activeSunspots.clear();
    for (int i = 0; i < count; ++i) {
        double theta = ((double)rand() / RAND_MAX) * 2.0 * 3.14159265358979;
        double phi = std::acos(2.0 * ((double)rand() / RAND_MAX) - 1.0);
        
        glm::dvec3 localX(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
        glm::dvec3 tempUp = (std::abs(localX.y) < 0.9) ? glm::dvec3(0.0, 1.0, 0.0) : glm::dvec3(0.0, 0.0, 1.0);
        glm::dvec3 localY = glm::normalize(glm::cross(tempUp, localX));
        glm::dvec3 localZ = glm::normalize(glm::cross(localX, localY));
        
        SunspotAnchor spot;
        spot.staticPos = localX * STAR_RADIUS;
        spot.toStaticTransform = glm::dmat3(localX, localY, localZ);
        spot.toLocalTransform = glm::transpose(spot.toStaticTransform);
        spot.tiltAngle = (((double)rand() / RAND_MAX) - 0.5) * 2.094; 

        spot.freeMagneticEnergy = 0.0;
        spot.trappedPlasmaMass = 1.0e10; // Base seed mass (kg)

        activeSunspots.push_back(spot);
    }
}

void updateCachedSunspotPos() {
    // Self-initialize random global prominence regions if none exist yet, like your-... nvm
    if (activeSunspots.empty()) {
        initializeRandomSunspots(5);
    }

    cachedCosRotation = std::cos(starRotationAngle);
    cachedSinRotation = std::sin(starRotationAngle);
    cachedCosInverseRotation = std::cos(-starRotationAngle);
    cachedSinInverseRotation = std::sin(-starRotationAngle);
    
    glm::dvec3 primaryStaticPos = activeSunspots[0].staticPos;
    cachedSunspotPos = glm::dvec3(
        primaryStaticPos.x * cachedCosRotation - primaryStaticPos.z * cachedSinRotation,
        primaryStaticPos.y,
        primaryStaticPos.x * cachedSinRotation + primaryStaticPos.z * cachedCosRotation
    );
}

glm::dvec3 calculateCurrentDensity(const glm::dvec3& pos){
    glm::dvec3 curlB = curl(pos, calculateMagneticField);
    return curlB / mu_0;
}

double calculateMagneticMoment(const glm::dvec3& velocity, double mass, const glm::dvec3& B){
    double B_mag = glm::length(B);
    if (B_mag < 1e-12) return 0.0;

    glm::dvec3 b_hat = B / B_mag;
    glm::dvec3 v_parallel = glm::dot(velocity, b_hat) * b_hat;
    glm::dvec3 v_perp = velocity - v_parallel;
    double v_perp_sq = glm::dot(v_perp, v_perp);

    double speed_sq = glm::dot(velocity, velocity);
    double gamma = 1.0 / std::sqrt(1.0 - speed_sq / (c * c));

    return (gamma * mass * v_perp_sq) / (2.0 * B_mag);
}

glm::dvec3 calculateGradBDrift(const glm::dvec3& pos, double charge, double mu) {
    glm::dvec3 B = calculateMagneticField(pos);
    double B_mag = glm::length(B);
    if (B_mag < 1e-12) return glm::dvec3(0.0);

    glm::dvec3 gradB = gradient(pos, getFieldMagnitude);
    return (mu / (charge * B_mag * B_mag)) * glm::cross(B, gradB);
}

double calculateMirrorForce(const glm::dvec3& pos, double mu) {
    glm::dvec3 B = calculateMagneticField(pos);
    double B_mag = glm::length(B);
    if (B_mag < 1e-12) return 0.0;

    glm::dvec3 b_hat = B / B_mag;
    glm::dvec3 gradB = gradient(pos, getFieldMagnitude);

    return -mu * glm::dot(b_hat, gradB);
}

double calculateAlfvenRadius() {
    static const double alfvenRadius = []() {
        double mu = globalStrength * std::pow(STAR_RADIUS, 3.0);
        return std::pow(mu, 4.0 / 7.0) 
             * std::pow(2.0 * G * STAR_MASS, -1.0 / 7.0) 
             * std::pow(MASS_ACCRETION_RATE, -2.0 / 7.0);
    }();
    return alfvenRadius;
}

glm::dvec3 calculateElectricField(const glm::dvec3& position, const glm::dvec3& B) {
    double r_len = glm::length(position);
    if (r_len < STAR_RADIUS) return glm::dvec3(0.0);

    // Heliocentric plasma framing, decoupling at the solar Alfven boundary
    double alfvenRadius = calculateAlfvenRadius(); 
    double decouplingFactor = std::exp(-(r_len - STAR_RADIUS) / alfvenRadius);

    glm::dvec3 rotationAxis = glm::dvec3(0.0, STAR_OMEGA, 0.0); 
    glm::dvec3 v_corotation = glm::cross(position, rotationAxis);

    glm::dvec3 radialDir = position / r_len;
    double latitudeFactor = std::abs(radialDir.y); 
    
    // Parker Hydrodynamic Acceleration Profile.... marshmellows
    double v_terminal = (0.0011 * c) + (latitudeFactor * (0.0013 * c));
    double scaleDistance = STAR_RADIUS * 4.0;
    double v_thermal = 0.00015 * c;

    double solarWindVel = v_thermal + (v_terminal * (1.0 - std::exp(-(r_len - STAR_RADIUS) / scaleDistance)));
    if (solarWindVel < 0.0) solarWindVel = 0.0;

    // Rotate query point backwards into the star's static frame to match B-field tracking 
    double cosA = cachedCosInverseRotation;
    double sinA = cachedSinInverseRotation;
    glm::dvec3 rotatedPos = glm::dvec3(
        position.x * cosA - position.z * sinA,
        position.y,
        position.x * sinA + position.z * cosA
    );

    // Accumulate the solar wind suppression factor across all localized anchors
    double combinedWindTaper = 1.0;
    double loopRadius = STAR_RADIUS * 0.65; 
    
    for (const auto& spot : activeSunspots) {
        double distToSunspot = glm::length(rotatedPos - spot.staticPos);
        if (distToSunspot < loopRadius) {
            double t = distToSunspot / loopRadius;
            double smoothTaper = 1.0 - (10.0 * t * t * t) + (15.0 * t * t * t * t) - (6.0 * t * t * t * t * t);
            combinedWindTaper *= (1.0 - smoothTaper); 
        }
    }
    solarWindVel *= combinedWindTaper;

    glm::dvec3 v_solarwind = radialDir * solarWindVel;

    // Ideal MHD Convective Field: E = -v x B (cross product)
    glm::dvec3 bulkPlasmaVelocity = (decouplingFactor * v_corotation) + ((1.0 - decouplingFactor) * v_solarwind);
    return -glm::cross(bulkPlasmaVelocity, B);
}

glm::dvec3 calculateMagneticField(const glm::dvec3& position) {
    //Rotate query point backwards into the star's static frame
    double cosA = cachedCosInverseRotation;  
    double sinA = cachedSinInverseRotation;
    
    glm::dvec3 rotatedPos = glm::dvec3(
        position.x * cosA - position.z * sinA,
        position.y,
        position.x * sinA + position.z * cosA
    );

    glm::dvec3 r = rotatedPos;    
    double r_len = glm::length(r);
    glm::dvec3 B_total = glm::dvec3(0.0);

    double tiltAngle = 0.1265; //in radians
    double alfvenRadius = calculateAlfvenRadius(); 
    
    if (r_len < STAR_RADIUS) {
        B_total = (2.0 * globalStrength / std::pow(STAR_RADIUS, 3.0)) * glm::dvec3(std::sin(tiltAngle), std::cos(tiltAngle), 0.0);
    } else {
        if (useParkerSpiral) {
            glm::dvec3 radialDir = r / r_len;
            double latitudeFactor = std::abs(radialDir.y); 
            
            double v_terminal = (0.0011 * c) + (latitudeFactor * (0.0013 * c));
            double scaleDistance = STAR_RADIUS * 4.0;
            double v_thermal = 0.00015 * c;
  
            double solarWindVel = v_thermal + (v_terminal * (1.0 - std::exp(-(r_len - STAR_RADIUS) / scaleDistance)));
            if (solarWindVel < 0.0) solarWindVel = 0.0;
            
            double transitLag = 0.0;
            if (r_len > alfvenRadius) {
                transitLag = STAR_OMEGA * (r_len - alfvenRadius) / solarWindVel;
            }
            
            double cosLag = std::cos(transitLag);
            double sinLag = std::sin(transitLag);
            glm::dvec3 laggedR = glm::dvec3(
                r.x * cosLag + r.z * sinLag,
                r.y,
                -r.x * sinLag + r.z * cosLag
            );

            double tiltedY = laggedR.x * std::sin(tiltAngle) + laggedR.y * std::cos(tiltAngle);
            double sheetThickness = STAR_RADIUS * 0.05; 
            double polarityBlend = std::tanh(tiltedY / sheetThickness);

            double B_r = globalStrength * (STAR_RADIUS * STAR_RADIUS) / (r_len * r_len) * polarityBlend;
            glm::dvec3 B_radial = radialDir * B_r; 
            
            glm::dvec3 rotationAxis = glm::dvec3(0.0, 1.0, 0.0);
            glm::dvec3 toroidalDir = glm::cross(radialDir, rotationAxis);
            double t_len = glm::length(toroidalDir);
            toroidalDir = (t_len > 0.001) ? (toroidalDir / t_len) : glm::dvec3(0.0);

            double B_phi = 0.0;
            if (r_len > alfvenRadius) {
                double sinTheta = glm::length(glm::dvec3(radialDir.x, 0.0, radialDir.z));

                double ratio = STAR_RADIUS / r_len;
                double correctionFactor = 1.0 - (ratio * ratio);
                
                B_phi = -B_r * ((STAR_OMEGA * r_len * sinTheta) / solarWindVel) * correctionFactor;
            }
            
            B_total = B_radial + (toroidalDir * B_phi);
        } else {
            glm::dvec3 m = glm::dvec3(std::sin(tiltAngle), std::cos(tiltAngle), 0.0) * globalStrength * std::pow(STAR_RADIUS, 3.0);
            B_total = (3.0 * glm::dot(m, r) * r - (r_len * r_len) * m) / std::pow(r_len, 5.0);
        }
    }

    // Localized Active Region 
    for (const auto& spot : activeSunspots) {
        glm::dvec3 localPos = spot.toLocalTransform * rotatedPos;
        glm::dvec3 localSunspotPos = glm::dvec3(STAR_RADIUS, 0.0, 0.0);
        
        // Dynamic loop radius swells upward with energy inflation
        double maxDist = STAR_RADIUS * 0.65 * spot.inflation; 
        glm::dvec3 deltaPos = localPos - localSunspotPos;
        double d = glm::length(deltaPos);
        
        if (d < maxDist) {
            double cosT = std::cos(spot.tiltAngle);
            double sinT = std::sin(spot.tiltAngle);

            // Inflate pole separation and height as energy builds up
            double poleDist = 0.18 * STAR_RADIUS * spot.inflation;

            glm::dvec3 poleOffset(0.0, -poleDist * sinT, poleDist * cosT);

            glm::dvec3 pos_N = localSunspotPos * 0.95 + poleOffset;
            glm::dvec3 pos_S = localSunspotPos * 0.95 - poleOffset;
            
            glm::dvec3 r_N = localPos - pos_N;
            glm::dvec3 r_S = localPos - pos_S;
            
            double len_N = std::max(glm::length(r_N), STAR_RADIUS * 0.01);
            double len_S = std::max(glm::length(r_S), STAR_RADIUS * 0.01);
            
            double spotStrength = globalStrength * 0.85 * std::pow(STAR_RADIUS, 2.0);

            // Standard Closed Dipole Loop
            glm::dvec3 B_closed = spotStrength * (r_N / std::pow(len_N, 3.0) - r_S / std::pow(len_S, 3.0));
            
            // Reconnected Open Radial Field Lines
            glm::dvec3 openDir = glm::normalize(localPos);
            glm::dvec3 B_open = openDir * (spotStrength / (len_N * len_N));

            // Magnetic Reconnection: Blend closed loops into open lines during eruption
            glm::dvec3 B_loop = glm::mix(B_closed, B_open, spot.eruptionFactor);
            
            double t = d / maxDist;
            double taper = 1.0 - (10.0 * t * t * t) + (15.0 * t * t * t * t) - (6.0 * t * t * t * t * t);
            
            B_total += (spot.toStaticTransform * B_loop) * taper;
        }
    }

    double cosF = cachedCosRotation;
    double sinF = cachedSinRotation;
    return glm::dvec3(
        B_total.x * cosF - B_total.z * sinF,
        B_total.y,
        B_total.x * sinF + B_total.z * cosF
    );
}

// Emitters & Injectors (Warm Thermal Solar Production)

void triggerCoronalMassEjection(std::vector<PlasmaParticle>& particles) {
    int particlesInBurst = 400; 
    particles.reserve(particles.size() + (particlesInBurst * 2)); 
    glm::dvec3 eruptionDir = glm::normalize(cachedSunspotPos); 

    for (int i = 0; i < particlesInBurst; ++i) {
        glm::dvec3 blastExpansion((double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5);
        glm::dvec3 bubbleDir = glm::normalize(eruptionDir * 1.2 + blastExpansion * 0.8);

        double cmeSpeed = (0.0033 * c) + ((double)rand() / RAND_MAX) * (0.0033 * c); 
        glm::dvec3 velocity = bubbleDir * cmeSpeed;
        glm::dvec3 spawnOffset = bubbleDir * (STAR_RADIUS * 0.02);

        PlasmaParticle electron;
        electron.position = cachedSunspotPos + spawnOffset;
        electron.velocity = velocity;
        electron.charge = Electron_Charge;
        electron.mass = ELECTRON_MASS;
        electron.color = glm::dvec3(0.5, 0.75, 0.99); 
        electron.lifetime = 12.0; 
        electron.active = true;
        electron.type = Particle_Type::ElectromagneticParticle;
        electron.escaped = false; 
        electron.gcaInitialized = false;
        particles.push_back(electron);

        PlasmaParticle ion;
        ion.position = cachedSunspotPos + spawnOffset;
        ion.velocity = velocity;
        ion.charge = Proton_Charge;
        ion.mass = PROTON_MASS;
        ion.color = glm::dvec3(1.0, 0.48, 0.08); 
        ion.lifetime = 12.0; 
        ion.active = true;
        ion.type = Particle_Type::ElectromagneticParticle;
        ion.escaped = false; 
        ion.gcaInitialized = false;
        particles.push_back(ion);
    }
}

void triggerCarringtonEvent(std::vector<PlasmaParticle>& particles) {
    int particlesInBurst = 1600; 
    particles.reserve(particles.size() + (particlesInBurst * 2)); 
    glm::dvec3 eruptionDir = glm::normalize(cachedSunspotPos);

    for (int i = 0; i < particlesInBurst; ++i) {
        glm::dvec3 blastExpansion((double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5);
        glm::dvec3 bubbleDir = glm::normalize(eruptionDir * 0.9 + blastExpansion * 1.0);

        double cmeSpeed = (0.008 * c) + ((double)rand() / RAND_MAX) * (0.002 * c);
        glm::dvec3 velocity = bubbleDir * cmeSpeed;
        glm::dvec3 spawnOffset = bubbleDir * (STAR_RADIUS * 0.04);

        PlasmaParticle electron;
        electron.position = cachedSunspotPos + spawnOffset;
        electron.velocity = velocity;
        electron.charge = Electron_Charge;
        electron.mass = ELECTRON_MASS;
        electron.color = glm::dvec3(1.0, 1.0, 0.85); 
        electron.lifetime = 8.0; 
        electron.active = true;
        electron.type = Particle_Type::ElectromagneticParticle;
        electron.escaped = false;
        electron.gcaInitialized = false;
        particles.push_back(electron);

        PlasmaParticle ion;
        ion.position = cachedSunspotPos + spawnOffset;
        ion.velocity = velocity;
        ion.charge = Proton_Charge;
        ion.mass = PROTON_MASS;
        ion.color = glm::dvec3(0.95, 0.28, 0.0); 
        ion.lifetime = 8.0; 
        ion.active = true;
        ion.type = Particle_Type::ElectromagneticParticle;
        ion.escaped = false;
        ion.gcaInitialized = false;
        particles.push_back(ion);
    }
}

void triggerLocalizedCME(std::vector<PlasmaParticle>& particles, const glm::dvec3& spotPos, double calculatedSpeed, int particleCount) {
    particles.reserve(particles.size() + (particleCount * 2)); 
    glm::dvec3 eruptionDir = glm::normalize(spotPos); 

    for (int i = 0; i < particleCount; ++i) {
        glm::dvec3 blastExpansion((double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5, (double)rand()/RAND_MAX - 0.5);
        glm::dvec3 bubbleDir = glm::normalize(eruptionDir * 2.0 + blastExpansion * 0.4);

        double speedVariation = calculatedSpeed * (0.85 + 0.3 * ((double)rand() / RAND_MAX)); 
        glm::dvec3 velocity = bubbleDir * speedVariation;

        // Small local variance around the prominence
        glm::dvec3 spawnOffset = blastExpansion * (STAR_RADIUS * 0.03);

        PlasmaParticle electron;
        electron.position = spotPos + spawnOffset;
        electron.velocity = velocity;
        electron.charge = Electron_Charge;
        electron.mass = ELECTRON_MASS;
        electron.color = glm::dvec3(0.5, 0.75, 0.99); 
        electron.lifetime = 12.0; 
        electron.active = true;
        electron.type = Particle_Type::ElectromagneticParticle;
        electron.escaped = false; 
        electron.gcaInitialized = false;
        particles.push_back(electron);

        PlasmaParticle ion;
        ion.position = spotPos + spawnOffset;
        ion.velocity = velocity;
        ion.charge = Proton_Charge;
        ion.mass = PROTON_MASS;
        ion.color = glm::dvec3(1.0, 0.48, 0.08); 
        ion.lifetime = 12.0; 
        ion.active = true;
        ion.type = Particle_Type::ElectromagneticParticle;
        ion.escaped = false; 
        ion.gcaInitialized = false;
        particles.push_back(ion);
    }
}

void injectPlasmaField(std::vector<PlasmaParticle>& particles) {
    const int particlesPerFrame = 20; 
    particles.reserve(particles.size() + (particlesPerFrame * 2)); 

    int ambientWindCount = static_cast<int>(particlesPerFrame * 0.75);
    int coronalLoopCount = particlesPerFrame - ambientWindCount;
    
    // Population A: Volumetric Coronal Ambient Wind
    for (int i = 0; i < ambientWindCount; ++i) {
        double theta = ((double)rand() / RAND_MAX) * 2.0 * 3.14159265359; //x,z plane
        double phi = std::acos(2.0 * ((double)rand() / RAND_MAX) - 1.0); //y height angle
        glm::dvec3 spawnDir(std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta), std::cos(phi));
        
        double latitudeFactor = std::abs(spawnDir.y); 
        double coronaExtension = STAR_RADIUS * (1.01 + ((double)rand() / RAND_MAX) * 0.12);
        glm::dvec3 basePosition = spawnDir * coronaExtension;

        double v_terminal = (0.0011 * c) + (latitudeFactor * (0.0013 * c));
        double scaleDistance = STAR_RADIUS * 4.0;
        double v_thermal = 0.00015 * c; 
        double windSpeed = v_thermal + (v_terminal * (1.0 - std::exp(-(coronaExtension - STAR_RADIUS) / scaleDistance)));
        if (windSpeed < 0.0) windSpeed = 0.0;

        glm::dvec3 radialNormal = spawnDir;
        glm::dvec3 tangent1 = glm::normalize(glm::cross(radialNormal, glm::dvec3(0.0, 1.0, 0.0)));
        glm::dvec3 tangent2 = glm::normalize(glm::cross(radialNormal, tangent1));

        // ~1-2 km/s characteristic supergranulation velocity in the transverse plane
        double turbulenceScale = 0.000005 * c; 
        double u_phi_fluctuation = (((double)rand() / RAND_MAX) - 0.5) * turbulenceScale;
        double u_theta_fluctuation = (((double)rand() / RAND_MAX) - 0.5) * turbulenceScale;

        glm::dvec3 turbulenceVector = (tangent1 * u_phi_fluctuation) + (tangent2 * u_theta_fluctuation);

        glm::dvec3 baseVelocity = (spawnDir * windSpeed) + turbulenceVector;

        glm::dvec3 e_offset(((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5));
        
        PlasmaParticle electron;
        electron.position = basePosition + e_offset * (STAR_RADIUS * 0.02);
        electron.velocity = baseVelocity + e_offset * (0.00005 * c);
        electron.charge = Electron_Charge;
        electron.mass = ELECTRON_MASS;
        electron.color = glm::dvec3(0.02, 0.12, 0.82); 
        electron.lifetime = 13.0; 
        electron.active = true;
        electron.type = Particle_Type::ElectromagneticParticle;
        electron.escaped = false;
        electron.gcaInitialized = false;
        particles.push_back(electron);

        glm::dvec3 i_offset(((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5));
        
        PlasmaParticle ion;
        ion.position = basePosition + i_offset * (STAR_RADIUS * 0.02);
        ion.velocity = baseVelocity + i_offset * (0.00005 * c);
        ion.charge = Proton_Charge;  
        ion.mass = PROTON_MASS;
        ion.color = glm::dvec3(0.84, 0.12, 0.02); 
        ion.lifetime = 13.0;
        ion.active = true;
        ion.type = Particle_Type::ElectromagneticParticle;
        ion.escaped = false;
        ion.gcaInitialized = false;
        particles.push_back(ion);

        /*glm::dvec3 he_offset(((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5), ((double)rand()/RAND_MAX - 0.5));

        PlasmaParticle He;
        He.position = basePosition + he_offset * (STAR_RADIUS * 0.02);
        He.velocity = baseVelocity + he_offset * (0.00005 * c);
        He.charge = Proton_Charge * (((double)rand() / RAND_MAX) > 0.5 ? 1.0 : 2.0);  
        He.mass = PROTON_MASS * 2.0;
        He.color = glm::dvec3(0.84, 0.12, 0.02); 
        He.lifetime = 16.0;
        He.active = true;
        He.type = Particle_Type::ElectromagneticParticle;
        He.escaped = false;
        ion.gcaInitialized = false;
        particles.push_back(He);*/
    }

    // Population B: Trapped Prominence / Active-Region Loops distributed across the anchor pool
    for (int i = 0; i < coronalLoopCount; ++i) {
        if (activeSunspots.empty()) break;
        
        // Cycle particles sequentially across active prominence regions
        const auto& spot = activeSunspots[i % activeSunspots.size()];
        
        double archAngle = ((double)rand() / RAND_MAX) * 3.14159265359;
        double archRadius = 0.18 * STAR_RADIUS; 
        double archHeight = 0.35 * STAR_RADIUS; 

        // Calculate offsets relative to local sunspot center
        double offsetRadial = 0.01 * STAR_RADIUS + std::sin(archAngle) * archHeight; 
        double offsetLateral = std::cos(archAngle) * archRadius; 

        double cosT = std::cos(spot.tiltAngle);
        double sinT = std::sin(spot.tiltAngle);

        glm::dvec3 localArchPos;
        localArchPos.x = STAR_RADIUS + offsetRadial;
        localArchPos.y = -offsetLateral * sinT;
        localArchPos.z =  offsetLateral * cosT;
        
        glm::dvec3 staticArchPos = spot.toStaticTransform * localArchPos;
        
        glm::dvec3 basePosition = glm::dvec3(
            staticArchPos.x * cachedCosRotation - staticArchPos.z * cachedSinRotation,
            staticArchPos.y,
            staticArchPos.x * cachedSinRotation + staticArchPos.z * cachedCosRotation
        );

        glm::dvec3 B_local = calculateMagneticField(basePosition);
        glm::dvec3 B_dir = glm::length(B_local) > 1e-10 ? glm::normalize(B_local) : glm::dvec3(0.0, 1.0, 0.0);
        
        // Decouple
        double e_flow = ((double)rand() / RAND_MAX > 0.5) ? 1.0 : -1.0;
        glm::dvec3 e_vel = B_dir * (0.0003 * c) * e_flow;
        e_vel += glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (0.00001 * c);

        PlasmaParticle electron;
        electron.position = basePosition + glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (STAR_RADIUS * 0.01);
        electron.velocity = e_vel;
        electron.charge = Electron_Charge;
        electron.mass = ELECTRON_MASS;
        electron.color = glm::dvec3(0.8, 0.5, 0.15); 
        electron.lifetime = 13.0; 
        electron.active = true;
        electron.type = Particle_Type::ElectromagneticParticle;
        electron.escaped = false;
        electron.gcaInitialized = false;
        particles.push_back(electron);

        // Decouple
        double i_flow = ((double)rand() / RAND_MAX > 0.5) ? 1.0 : -1.0;
        glm::dvec3 i_vel = B_dir * (0.0003 * c) * i_flow;
        i_vel += glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (0.00001 * c);

        PlasmaParticle ion;
        ion.position = basePosition + glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (STAR_RADIUS * 0.01);
        ion.velocity = i_vel;
        ion.charge = Proton_Charge;
        ion.mass = PROTON_MASS;
        ion.color = glm::dvec3(0.95, 0.7, 0.02); 
        ion.lifetime = 13.0;
        ion.active = true;
        ion.type = Particle_Type::ElectromagneticParticle;
        ion.escaped = false;
        ion.gcaInitialized = false;
        particles.push_back(ion);

        double He_flow = ((double)rand() / RAND_MAX > 0.5) ? 1.0 : -1.0;
        glm::dvec3 He_vel = B_dir * (0.0003 *c ) * He_flow;
        He_vel += glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (0.00001 * c);

        PlasmaParticle He;// would either be He3 or He4, with charges 1 -> 2
        He.position = basePosition + glm::dvec3((double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5, (double)rand()/RAND_MAX-0.5) * (STAR_RADIUS * 0.01);
        He.velocity = He_vel;
        double HeState = ((double)rand() / RAND_MAX> 0.5) ? 1.0 : 2.0;
        He.charge =  Proton_Charge * HeState;
        He.mass = PROTON_MASS * 2 + NEUTRON_MASS * (1 + (rand() % 2));
        He.color = glm::dvec3(1.0, 0.71, 0.75);
        He.lifetime = 13.0;
        He.active = true;
        He.type = Particle_Type::ElectromagneticParticle;
        He.escaped = false;
        He.gcaInitialized = false;
        particles.push_back(He);

    }
}

void injectCosmicRay(std::vector<PlasmaParticle>& particles) {
    const int raysPerTrigger = 6;
    particles.reserve(particles.size() + raysPerTrigger); 
    double spawnRadius = SIM_LIMIT - (STAR_RADIUS * 1.5); 

    for (int i = 0; i < raysPerTrigger; ++i) {
        PlasmaParticle ray;
        double theta = ((double)rand() / RAND_MAX) * 2.0 * 3.14159265359;
        double phi = std::acos(2.0 * ((double)rand() / RAND_MAX) - 1.0);

        ray.position.x = spawnRadius * std::sin(phi) * std::cos(theta);
        ray.position.y = spawnRadius * std::sin(phi) * std::sin(theta);
        ray.position.z = spawnRadius * std::cos(phi);

        glm::dvec3 target = glm::dvec3(((double)rand()/RAND_MAX - 0.5)*STAR_RADIUS, ((double)rand()/RAND_MAX - 0.5)*STAR_RADIUS, ((double)rand()/RAND_MAX - 0.5)*STAR_RADIUS);
        double energyRoll = (double)rand() / RAND_MAX;
        bool isUHECR = (energyRoll < 0.02); 

        if (isUHECR) {
            ray.velocity = glm::normalize(target - ray.position) * ((0.99 * c) + ((double)rand() / RAND_MAX) * (0.009 * c));
            ray.lifetime = 1.5;
        } else {
            ray.velocity = glm::normalize(target - ray.position) * ((0.80 * c) + ((double)rand() / RAND_MAX) * (0.10 * c));
            ray.lifetime = 4.5;
        }
        
        double chargeRoll = (double)rand() / RAND_MAX;
        if (chargeRoll < 0.89) {
            ray.charge = Proton_Charge; ray.mass = PROTON_MASS;
            ray.color = isUHECR ? glm::dvec3(1.0, 0.0, 1.0) : glm::dvec3(0.3, 1.0, 0.3);
        } else if (chargeRoll < 0.98) {
            ray.charge = Proton_Charge * 2.0; ray.mass = PROTON_MASS * 4.0;
            ray.color = isUHECR ? glm::dvec3(1.0, 0.2, 1.0) : glm::dvec3(0.8, 0.7, 0.0);
        } else {
            ray.charge = Electron_Charge; ray.mass = ELECTRON_MASS;
            ray.color = glm::dvec3(0.5, 0.6, 1.0);
        }

        ray.active = true;
        ray.type = Particle_Type::CosmicRayParticle;
        ray.escaped = true;
        ray.gcaInitialized = false;
        particles.push_back(ray);
    }
}

void handleCollisions(std::vector<PlasmaParticle>& particles) {

} // place holder,,, idk but this might be useful later on

int identifySourceSunspot(const PlasmaParticle& p) {
    double r_len = glm::length(p.position);
    if (r_len < STAR_RADIUS * 2.0) return -1; 

    // Convert Cartesian to Spherical (Helio-equatorial)
    double theta = std::asin(p.position.y / r_len); // Latitude
    double phi = std::atan2(p.position.z, p.position.x); // Longitude
    
    // Estimate average solar wind speed for this latitude
    double latitudeFactor = std::abs(p.position.y / r_len);
    double v_terminal = (0.0011 * c) + (latitudeFactor * (0.0013 * c));
    double v_sw = v_terminal; // Approximate bulk average transit speed

    double transitTime = (r_len - STAR_RADIUS) / v_sw;
    double phaseLag = STAR_OMEGA * transitTime;
    
    // Theoretical Source Longitude
    double sourcePhi = phi + phaseLag; 
    
    // Convert source spherical coordinates back to star's local Cartesian frame
    double cosLat = std::cos(theta);
    glm::dvec3 theoreticalSource(
        STAR_RADIUS * cosLat * std::cos(sourcePhi),
        STAR_RADIUS * std::sin(theta),
        STAR_RADIUS * cosLat * std::sin(sourcePhi)
    );

    // Rotate theoretical source into the static anchor frame based on current star Rotation Angle
    double cosF = std::cos(-starRotationAngle);
    double sinF = std::sin(-starRotationAngle);
    glm::dvec3 staticSource(
        theoreticalSource.x * cosF - theoreticalSource.z * sinF,
        theoreticalSource.y,
        theoreticalSource.x * sinF + theoreticalSource.z * cosF
    );

    // Finds closest active sunspot anchor
    int bestMatch = -1;
    double minDistance = std::numeric_limits<double>::max();
    
    for (int i = 0; i < activeSunspots.size(); ++i) {
        double dist = glm::length(activeSunspots[i].staticPos - staticSource);
        if (dist < minDistance) {
            minDistance = dist;
            bestMatch = i;
        }
    }
    
    return bestMatch;
}

void evaluateMHDInstabilities(std::vector<PlasmaParticle>& particles, double physicsDt) {
    if (activeSunspots.empty()) return;

    double loopRadius = STAR_RADIUS * 0.35;
    double loopVolume = (4.0 / 3.0) * 3.14159265359 * std::pow(loopRadius, 3.0);
    double B_spot = globalStrength * 0.85; 
    
    double E_potential = (B_spot * B_spot / (2.0 * mu_0)) * loopVolume;

    for (auto& spot : activeSunspots) {
        double latitudeFactor = std::abs(spot.staticPos.y / STAR_RADIUS);
        
        double shearRate = 2.5e-1 * STAR_OMEGA * (1.0 + latitudeFactor); 
        spot.freeMagneticEnergy += E_potential * shearRate * physicsDt;

        double massAccretionRate = 1.0e8;
        spot.trappedPlasmaMass += massAccretionRate * physicsDt;

        double eta = spot.freeMagneticEnergy / E_potential;

        //loop swells up to 80% higher prior to collapse
        spot.inflation = 1.0 + 0.8 * std::min(eta / 0.28, 1.0);

        // decay open field back to closed loops over time
        if (spot.eruptionFactor > 0.0) {
            spot.eruptionFactor = std::max(0.0, spot.eruptionFactor - 0.2 * physicsDt);
        }

        double randomEtaThreshold = 0.28 + (((double)rand() / RAND_MAX) - 0.5) * 0.05;

        if (eta >= randomEtaThreshold) {
            double v_cme = std::sqrt((2.0 * spot.freeMagneticEnergy) / spot.trappedPlasmaMass);
            double speedVariance = 0.85 + 0.3 * ((double)rand() / RAND_MAX);
            v_cme = std::clamp(v_cme * speedVariance, 100.0e3, 3000.0e3);

            // Trigger Magnetic Reconnection
            spot.eruptionFactor = 1.0; 

            // Spawn CME from the actual inflated prominence height
            glm::dvec3 currentSpotPos = glm::dvec3(
                spot.staticPos.x * cachedCosRotation - spot.staticPos.z * cachedSinRotation,
                spot.staticPos.y,
                spot.staticPos.x * cachedSinRotation + spot.staticPos.z * cachedCosRotation
            );
            glm::dvec3 prominencePos = currentSpotPos * (1.0 + 0.15 * spot.inflation);

            int particleBurstCount = static_cast<int>(300 * (eta / 0.28));

            std::cout << "[MHD RECONNECTION] Prominence collapsed at latitude factor: " 
                      << latitudeFactor << " | Speed: " << (v_cme / 1000.0) << " km/s" << std::endl;

            triggerLocalizedCME(particles, prominencePos, v_cme, particleBurstCount);

            // Reset energy state
            spot.freeMagneticEnergy = 0.0;
            spot.trappedPlasmaMass = 1.0e10; 
            spot.inflation = 1.0; 
        }
    }
}

// Interaction / Integration Functions 

void updateChargedParticle(PlasmaParticle& p, double dt, double starRadius) {
    double velocity_squared = glm::dot(p.velocity, p.velocity);
    double c_squared = c * c; 
    if (velocity_squared >= c_squared * 0.9999) velocity_squared = c_squared * 0.9999;

    double effectiveMass = p.mass; 
    double chargeOverMass = p.charge / effectiveMass;

    glm::dvec3 Magnetic_field = calculateMagneticField(p.position);
    double B_mag = glm::length(Magnetic_field);
    
    double gamma = 1.0 / std::sqrt(1.0 - (velocity_squared / c_squared));
    double omega_c = (std::abs(p.charge) * B_mag) / (effectiveMass * gamma);
    
    if (omega_c * dt > 25.0) {
        glm::dvec3 Electric_field = enableElectricField ? calculateElectricField(p.position, Magnetic_field) : glm::dvec3(0.0);
        double B_mag_sq = B_mag * B_mag;
        
        if (B_mag_sq > 1e-20) {
            glm::dvec3 b_hat = Magnetic_field / B_mag;
            glm::dvec3 v_drift = glm::cross(Electric_field, Magnetic_field) / B_mag_sq;

            if (!p.gcaInitialized) {
                glm::dvec3 v_perp = p.velocity - (b_hat * glm::dot(p.velocity, b_hat)) - v_drift;
                double v_perp_sq = glm::dot(v_perp, v_perp);
                double gamma = 1.0 / std::sqrt(1.0 - (velocity_squared / (c * c)));
                p.mu = (gamma * effectiveMass * v_perp_sq) / (2.0 * B_mag);
                p.v_parallel = glm::dot(p.velocity, b_hat);
                p.gcaInitialized = true;
            }

            double ds = starRadius * 0.002; 
            double B_plus = glm::length(calculateMagneticField(p.position + b_hat * ds));
            double B_minus = glm::length(calculateMagneticField(p.position - b_hat * ds));
            double grad_parallel_B = (B_plus - B_minus) / (2.0 * ds);

            glm::dvec3 grad_B(0.0);
            for(int i = 0; i < 3; ++i) {
                glm::dvec3 offset(0.0); offset[i] = ds;
                double B_p = glm::length(calculateMagneticField(p.position + offset));
                double B_m = glm::length(calculateMagneticField(p.position - offset));
                grad_B[i] = (B_p - B_m) / (2.0 * ds);
            }

            double acc_parallel = -(p.mu / effectiveMass) * grad_parallel_B;
            p.v_parallel += acc_parallel * dt;

            glm::dvec3 v_grad_drift = (p.mu / (p.charge * B_mag)) * glm::cross(b_hat, grad_B);

            p.velocity = (b_hat * p.v_parallel) + v_drift + v_grad_drift;
            
            double v_len_sq = glm::dot(p.velocity, p.velocity);
            if (v_len_sq >= c_squared) p.velocity = glm::normalize(p.velocity) * (0.999 * c);
            
            p.position += p.velocity * dt;

            double r_len_squared = glm::dot(p.position, p.position);
            if (r_len_squared >= (starRadius * starRadius)) p.escaped = true;
            if (r_len_squared < (starRadius * starRadius * 0.995)) p.active = false;
            return;
        }
    } else {
        p.gcaInitialized = false;
    }

    int subSteps = static_cast<int>(std::ceil(omega_c * dt / 0.5));
    subSteps = std::clamp(subSteps, 1, 32); 

    double subDt = dt / static_cast<double>(subSteps);
    double halfSubDt = subDt * 0.5;

    for (int step = 0; step < subSteps; ++step) {
        Magnetic_field = calculateMagneticField(p.position);
        glm::dvec3 Electric_field = enableElectricField ? calculateElectricField(p.position, Magnetic_field) : glm::dvec3(0.0);
        
        double velocity_squared_loop = glm::dot(p.velocity, p.velocity);
        if (velocity_squared_loop >= c_squared * 0.9999) velocity_squared_loop = c_squared * 0.9999;

        double gamma_loop = 1.0 / std::sqrt(1.0 - (velocity_squared_loop / c_squared));
        glm::dvec3 u = gamma_loop * p.velocity;

        glm::dvec3 E_accel = (chargeOverMass * Electric_field) * halfSubDt;
        glm::dvec3 u_minus = u + E_accel;

        double u_minus_sq = glm::dot(u_minus, u_minus);
        double gamma_minus = std::sqrt(1.0 + u_minus_sq / c_squared);

        glm::dvec3 t = (chargeOverMass / gamma_minus) * Magnetic_field * halfSubDt;
        glm::dvec3 s = (2.0 * t) / (1.0 + glm::dot(t, t));
        glm::dvec3 u_prime = u_minus + glm::cross(u_minus, t);
        glm::dvec3 u_plus = u_minus + glm::cross(u_prime, s);
        
        glm::dvec3 u_next = u_plus + E_accel;
        double gamma_next = std::sqrt(1.0 + glm::dot(u_next, u_next) / c_squared);
        
        p.velocity = u_next / gamma_next;
        p.position += p.velocity * subDt;

        double r_len_squared = glm::dot(p.position, p.position);
        if (r_len_squared >= (starRadius * starRadius)) p.escaped = true;
        if (r_len_squared < (starRadius * starRadius * 0.995)) {
            p.active = false;
            break;
        }
    }
}

void updateParticles(std::vector<PlasmaParticle>& particles, double dt) {
    const double simLimitSq = SIM_LIMIT * SIM_LIMIT;
    double timeScale = 6900.0;
    double physicsDt = dt * timeScale;

    starRotationAngle += STAR_OMEGA * physicsDt; 
    updateCachedSunspotPos(); 

    evaluateMHDInstabilities(particles, physicsDt);

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        auto& p = particles[i];
        if (!p.active) continue;

        p.lifetime -= dt;
        if (p.lifetime <= 0.0) { p.active = false; continue; }

        double r_len_squared = glm::dot(p.position, p.position);
        if (r_len_squared > simLimitSq) { p.active = false; continue; }

        if (p.charge != 0.0) { 
            updateChargedParticle(p, physicsDt, STAR_RADIUS);
        } else {
            p.position += p.velocity * physicsDt;
        }
    }

    // Fast swap-and-pop removal: O(N) linear scan without memory shifts
    size_t i = 0;
    while (i < particles.size()) {
        if (!particles[i].active) {
            particles[i] = particles.back();
            particles.pop_back();
        } else {
            ++i;
        }
    }
}