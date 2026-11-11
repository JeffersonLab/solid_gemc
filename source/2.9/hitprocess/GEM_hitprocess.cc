// Small-prototype Hall C GEM hit process.
//
// Adapted from the solid_gemc uRwell hit-process structure, but the
// large-uRwell stereo-strip assumption is removed.
//
// GEM readout:
//   layer 1 -> direct local X -> component 1...256
//   layer 2 -> direct local Y -> component 1...256
//
// No strip-angle rotation is used.
//
// Detector-response constants below are intentionally kept at the current
// uRwell values until GEM-specific values are confirmed.

#include "GEM_hitprocess.h"
#include "GEM_strip.h"

#include "G4Trap.hh"
#include "G4Box.hh"

#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

#include <iostream>
using namespace std;

static GEMConstants initializeGEMConstants(
    int runno,
    string digiVariation = "default",
    string digiSnapshotTime = "no",
    bool accountForHardwareStatus = false)
{
    GEMConstants c;

    if(runno == -1) return c;

    c.runNo = runno;
    c.date  = "2026-10-26";

    if(getenv("CCDB_CONNECTION") != nullptr)
        c.connection = (string)getenv("CCDB_CONNECTION");
    else
        c.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

    // ------------------------------------------------------------
    // Small GEM prototype: direct X-Y readout, 256 strips each.
    // ------------------------------------------------------------
    c.number_of_x_strips = 256;
    c.number_of_y_strips = 256;

    // Dimensions are filled from the actual Geant4 sensitive volume
    // in processID before FindStrip is called.
    c.Xhalf = 0.0;
    c.Yhalf = 0.0;
    c.Zhalf = 0.0;

    // ------------------------------------------------------------
    // TEMPORARY detector-response parameters inherited from uRwell.
    // Confirm/update GEM-specific values separately.
    // ------------------------------------------------------------
    c.w_i        = 25.0;  // eV
    c.sigma_td   = 0.5;
    c.nb_sigma   = 5;
    c.gain       = 1E4;
    c.v_drift    = 5E-3;  // cm/ns
    c.sigma_time = 20.0;  // ns

    return c;
}

map<string, double> GEM_HitProcess::integrateDgt(MHit* aHit, int hitn)
{
    map<string, double> dgtz;
    vector<identifier> identity = aHit->GetId();

    rejectHitConditions = false;
    writeHit = true;

    trueInfos tInfos(aHit);

    dgtz["hitn"]      = hitn;
    dgtz["sector"]    = identity[1].id;
    dgtz["layer"]     = identity[3].id;
    dgtz["component"] = identity[4].id;

    if(identity[4].id == -15000) {
        dgtz["ADC_ADC"]  = 0;
        dgtz["ADC_time"] = 0;
    }
    else {
        // Keep the current uRwell ADC convention.
        dgtz["ADC_ADC"] =
            (1.0 * (int)(gemC.gain * 1e6 * tInfos.eTot / gemC.w_i));
        dgtz["ADC_time"] = identity[4].time;
    }

    dgtz["ADC_ped"] = 0;

    if(rejectHitConditions)
        writeHit = false;

    return dgtz;
}

vector<identifier> GEM_HitProcess::processID(
    vector<identifier> id,
    G4Step* aStep,
    detector Detector)
{
    vector<identifier> yid;
    GEM_strip gemStrip;

    // Global -> local sensitive-volume position.
    G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();

    G4ThreeVector lxyz =
        aStep->GetPreStepPoint()
             ->GetTouchableHandle()
             ->GetHistory()
             ->GetTopTransform()
             .TransformPoint(xyz);

    G4VTouchable* TH =
        (G4VTouchable*)aStep->GetPreStepPoint()->GetTouchable();

    G4VSolid* solid = TH->GetSolid();

    G4Trap* Trap = dynamic_cast<G4Trap*>(solid);
    G4Box*  Box  = dynamic_cast<G4Box*>(solid);

    // Small Hall C GEMs are expected to use G4Box.
    // G4Trap support is retained to stay close to the original framework.
    if(Trap != nullptr) {
        gemC.Xhalf = Trap->GetXHalfLength1();
        gemC.Yhalf = Trap->GetYHalfLength1();
        gemC.Zhalf = Trap->GetZHalfLength();
    }
    else if(Box != nullptr) {
        gemC.Xhalf = Box->GetXHalfLength();
        gemC.Yhalf = Box->GetYHalfLength();
        gemC.Zhalf = Box->GetZHalfLength();
    }
    else {
        cout << "Warning: GEM sensitive volume is neither G4Trap nor G4Box."
             << endl;

        if(!id.empty())
            id[id.size()-1].id_sharing = 1.0;

        return id;
    }

    const double depe = aStep->GetTotalEnergyDeposit();
    const double time = aStep->GetPostStepPoint()->GetGlobalTime();

    if(depe <= 0.0) {
        if(!id.empty())
            id[id.size()-1].id_sharing = 1.0;

        return id;
    }

    // ============================================================
    // X READOUT
    // Direct local x coordinate; NO strip-angle rotation.
    // layer = 1
    // component = 1 ... 256
    // ============================================================
    vector<GEM_strip_found> hits_x =
        gemStrip.FindStrip(lxyz, depe, gemC, time, 'x');

    for(unsigned int h = 0; h < hits_x.size(); h++) {
        if(hits_x.at(h).numberID < 1) continue;

        for(unsigned int j = 0; j < id.size(); j++) {
            identifier this_id;

            this_id.name = id[j].name;
            this_id.rule = id[j].rule;
            this_id.id   = id[j].id;
            this_id.time = id[j].time;

            if(this_id.name == "layer")
                this_id.id = 1;

            if(this_id.name == "component") {
                this_id.id   = hits_x.at(h).numberID;
                this_id.time = hits_x.at(h).time;
            }

            this_id.TimeWindow = id[j].TimeWindow;
            this_id.TrackId    = id[j].TrackId;

            this_id.id_sharing =
                hits_x.at(h).weight /
                (gemC.gain * 1e6 * depe / gemC.w_i);

            yid.push_back(this_id);
        }
    }

    // ============================================================
    // Y READOUT
    // Direct local y coordinate; NO strip-angle rotation.
    // layer = 2
    // component = 1 ... 256
    // ============================================================
    vector<GEM_strip_found> hits_y =
        gemStrip.FindStrip(lxyz, depe, gemC, time, 'y');

    for(unsigned int h = 0; h < hits_y.size(); h++) {
        if(hits_y.at(h).numberID < 1) continue;

        for(unsigned int j = 0; j < id.size(); j++) {
            identifier this_id;

            this_id.name = id[j].name;
            this_id.rule = id[j].rule;
            this_id.id   = id[j].id;
            this_id.time = id[j].time;

            if(this_id.name == "layer")
                this_id.id = 2;

            if(this_id.name == "component") {
                this_id.id   = hits_y.at(h).numberID;
                this_id.time = hits_y.at(h).time;
            }

            this_id.TimeWindow = id[j].TimeWindow;
            this_id.TrackId    = id[j].TrackId;

            this_id.id_sharing =
                hits_y.at(h).weight /
                (gemC.gain * 1e6 * depe / gemC.w_i);

            yid.push_back(this_id);
        }
    }

    return yid;
}

vector<MHit*> GEM_HitProcess::electronicNoise()
{
    vector<MHit*> noiseHits;
    return noiseHits;
}

map<int, vector<double> > GEM_HitProcess::chargeTime(MHit* aHit, int hitn)
{
    map<int, vector<double> > CT;
    return CT;
}

double GEM_HitProcess::voltage(double charge, double time, double forTime)
{
    return 0.0;
}

map<string, vector<int> > GEM_HitProcess::multiDgt(MHit* aHit, int hitn)
{
    map<string, vector<int> > MH;
    return MH;
}

void GEM_HitProcess::initWithRunNumber(int runno)
{
    string digiVariation =
        gemcOpt.optMap["DIGITIZATION_VARIATION"].args;

    string digiSnapshotTime =
        gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;

    if(gemC.runNo != runno) {
        cout << " > Initializing " << HCname
             << " digitization for run number "
             << runno << endl;

        gemC = initializeGEMConstants(
            runno,
            digiVariation,
            digiSnapshotTime,
            accountForHardwareStatus);

        gemC.runNo = runno;
    }
}

GEMConstants GEM_HitProcess::gemC = initializeGEMConstants(1);
