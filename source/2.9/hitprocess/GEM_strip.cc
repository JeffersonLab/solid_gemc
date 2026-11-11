// Small-prototype GEM strip finding.
// Direct X-Y readout: no stereo-strip rotation.

#include "GEM_strip.h"

#include "Randomize.hh"
#include "G4Poisson.hh"

#include <cmath>

vector<GEM_strip_found> GEM_strip::FindStrip(G4ThreeVector xyz,
                                             double Edep,
                                             const GEMConstants& gemC,
                                             double time,
                                             char axis)
{
    vector<GEM_strip_found> strip_found;
    GEM_strip_found thisStrip;

    // Convert deposited energy (MeV) to primary electrons.
    int N_el = static_cast<int>(1e6 * Edep / gemC.w_i);

    if(N_el <= 0) {
        thisStrip.numberID = -15000;
        thisStrip.weight   = 1.0;
        thisStrip.time     = -1.0;
        strip_found.push_back(thisStrip);
        return strip_found;
    }

    // Keep the same gain treatment as the current uRwell implementation.
    N_el = G4Poisson(N_el * gemC.gain);

    double coordinate = 0.0;
    double halfSize   = 0.0;
    int nStrips       = 0;

    if(axis == 'x' || axis == 'X') {
        coordinate = xyz.x();
        halfSize   = gemC.Xhalf;
        nStrips    = gemC.number_of_x_strips;
    }
    else if(axis == 'y' || axis == 'Y') {
        coordinate = xyz.y();
        halfSize   = gemC.Yhalf;
        nStrips    = gemC.number_of_y_strips;
    }
    else {
        thisStrip.numberID = -15000;
        thisStrip.weight   = 1.0;
        thisStrip.time     = -1.0;
        strip_found.push_back(thisStrip);
        return strip_found;
    }

    // Direct X-Y strip mapping.
    // For 256 strips:
    //   coordinate = -halfSize  -> component 1
    //   coordinate ~ 0         -> component ~128/129
    //   coordinate < +halfSize -> component 256
    const double pitch = 2.0 * halfSize / nStrips;
    const int stripID =
        static_cast<int>(floor((coordinate + halfSize) / pitch)) + 1;

    if(stripID < 1 || stripID > nStrips) {
        thisStrip.numberID = -15000;
        thisStrip.weight   = 1.0;
        thisStrip.time     = -1.0;
        strip_found.push_back(thisStrip);
        return strip_found;
    }

    // Same drift-time prescription as the current uRwell code.
    const double z_real = xyz.z();
    const double time_dz =
        fabs(-gemC.Zhalf/cm + z_real/cm) / gemC.v_drift;

    thisStrip.numberID = stripID;
    thisStrip.weight   = N_el;
    thisStrip.time =
        time + time_dz + G4RandGauss::shoot(0.0, gemC.sigma_time);

    strip_found.push_back(thisStrip);
    return strip_found;
}
