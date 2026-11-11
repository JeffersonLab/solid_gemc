#ifndef GEM_STRIP_H
#define GEM_STRIP_H 1

#include "G4ThreeVector.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <string>
#include <vector>

using namespace CLHEP;
using namespace std;

// Constants used by the small-prototype GEM strip response.
// Detector-response values are presently inherited from the uRwell code
// and should be updated if GEM-specific values are established.
class GEMConstants
{
public:
    int runNo;
    string date;
    string connection;
    char database[80];

    // Sensitive-volume half dimensions, filled from Geant4 in processID.
    double Xhalf;
    double Yhalf;
    double Zhalf;

    // Small prototype GEM: 256 X strips and 256 Y strips.
    int number_of_x_strips;
    int number_of_y_strips;

    // Detector-response parameters.
    double w_i;
    double sigma_td;
    int nb_sigma;
    double gain;
    double v_drift;
    double sigma_time;
    double v_eff_readout = 0;
};

class GEM_strip_found
{
public:
    int numberID;
    double weight;
    double time;
};

class GEM_strip
{
public:
    // axis = 'x': direct local-x -> X strip
    // axis = 'y': direct local-y -> Y strip
    // No stereo-angle rotation is applied.
    vector<GEM_strip_found> FindStrip(G4ThreeVector xyz,
                                      double Edep,
                                      const GEMConstants& gemC,
                                      double time,
                                      char axis);
};

#endif
