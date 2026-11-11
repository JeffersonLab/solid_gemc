#ifndef GEM_HITPROCESS_H
#define GEM_HITPROCESS_H 1

#include "HitProcess.h"
#include "GEM_strip.h"

// Hit process for the small Hall C GEM prototypes.
// Readout is direct X-Y (0-degree, no stereo-strip rotation),
// with 256 strips in each coordinate.
class GEM_HitProcess : public HitProcess
{
public:
    ~GEM_HitProcess() {;}

    map<string, double> integrateDgt(MHit*, int);
    map<string, vector<int> > multiDgt(MHit*, int);
    virtual map<int, vector<double> > chargeTime(MHit*, int);
    virtual double voltage(double, double, double);

    vector<identifier> processID(vector<identifier>, G4Step*, detector);

    static HitProcess *createHitClass() { return new GEM_HitProcess; }

private:
    static GEMConstants gemC;
    void initWithRunNumber(int runno);
    vector<MHit*> electronicNoise();
};

#endif
