// gemc headers
#include "uRwell_hitprocess.h"
#include "uRwell_strip.h"

// geant4 headers
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4CachedMagneticField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "G4Trap.hh"
#include "G4Box.hh"


// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static uRwellConstants initializeuRwellConstants(int runno, string digiVariation = "default", string digiSnapshotTime = "no", bool accountForHardwareStatus = false)
{
	// all these constants should be read from CCDB
	uRwellConstants urwellC;
	
	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return urwellC;
	string timestamp = "";
	if(digiSnapshotTime != "no") {
		timestamp = ":"+digiSnapshotTime;
	}
	
	// database
	urwellC.runNo = runno;
	urwellC.date       = "2022-08-23";
	if(getenv ("CCDB_CONNECTION") != nullptr)
		urwellC.connection = (string) getenv("CCDB_CONNECTION");
	else
		urwellC.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";
	
	/*
	 unique_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(fmtc.connection));
	 vector<vector<double> > data;
	 //Load the geometrical constant for all layers
	 sprintf(fmtc.database,"/geometry/uRwell/uRwell_global:%d:%s%s", fmtc.runNo, digiVariation.c_str(), timestamp.c_str());
	 data.clear(); calib->GetCalib(data,uRwell.database);
	 // all dimensions are in mm
	 */
	/*number of strip in each chambers*/
	urwellC.number_strip_chamber[0] = 128;
	urwellC.number_strip_chamber[1] = 128;
	urwellC.number_strip_chamber[2] = 128;
	
	urwellC.number_of_strip = 128; //Total number of strip
	urwellC.stripU_stereo_angle = -10 ; // angle between strip and trapezoid base in degree
	urwellC.stripU_pitch = 1.;  //mm

	urwellC.stripU_width[0] = 0.4;
	urwellC.stripU_width[1] = 0.4;
	urwellC.stripU_width[2] = 0.4;
	urwellC.stripU_width[3] = 0.4;

	urwellC.stripU_width_proto[0] = 0.175;
	urwellC.stripU_width_proto[1] = 0.350;
	urwellC.stripU_width_proto[2] = 0.262;
	urwellC.stripU_width_proto[3] = 0.350;

	urwellC.stripV_stereo_angle = 10 ; // angle between strip and trapezoid base in degree
	urwellC.stripV_pitch = 1.;  //mm

	urwellC.stripV_width[0] = 0.4;
	urwellC.stripV_width[1] = 0.4;
	urwellC.stripV_width[2] = 0.4;
	urwellC.stripV_width[3] = 0.4;

	urwellC.stripV_width_proto[0] = 0.355;
	urwellC.stripV_width_proto[1] = 0.650;
	urwellC.stripV_width_proto[2] = 0.5;
	urwellC.stripV_width_proto[3] = 0.650;
	
	urwellC.w_i=25; //ionization potential assumed to be 25 eV
	urwellC.sigma_td= 0.5;         // effective value to take into account transverse diffusion + charge dispersion
	urwellC.nb_sigma = 5;            // Number of sigma to study around the closest strip
	urwellC.gain =1E4;
	
	// drift velocity
	urwellC.v_drift = 5E-3; // velocity drift [cm/ns]
	urwellC.sigma_time = 20; // time resolution 20 ns
	
	return urwellC;
}




map<string, double>uRwell_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	rejectHitConditions = false;
	writeHit = true;

	//uRwellConstants uRwellC;
	trueInfos tInfos(aHit);
        dgtz["hitn"]      = hitn;
        dgtz["sector"]    = identity[1].id;
        dgtz["layer"]     = identity[3].id;
        dgtz["component"] = identity[4].id;	
	
    if(identity[4].id ==-15000){
    	dgtz["ADC_ADC"]  = 0;
    	dgtz["ADC_time"] = 0;
    }else{
    	dgtz["ADC_ADC"]  = (1.0*(int) (uRwellC.gain*1e6*tInfos.eTot/uRwellC.w_i));
    	dgtz["ADC_time"] = identity[4].time;
    }
	dgtz["ADC_ped"]   = 0;
	
//	cout<<dgtz["sector"]<<" "<<dgtz["layer"]<<" "<<dgtz["component"]<<" "<<dgtz["ADC_ADC"]<<endl;
	
	// define conditions to reject hit
	if (rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
	
}

vector<identifier> uRwell_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
    vector<identifier> yid;
    uRwell_strip URwell_strip;

    // Local coordinates of the hit in the sensitive gas volume
    G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector lxyz =
        aStep->GetPreStepPoint()
             ->GetTouchableHandle()
             ->GetHistory()
             ->GetTopTransform()
             .TransformPoint(xyz);

    // The original LMU strip algorithm was written for G4Trap volumes.
    // For the HallC2026 compact rectangular geometries, the sensitive gas
    // volume can be a G4Box.  Treat a Box as a special trapezoid with
    // equal small and large x half-lengths, so the same strip-finding code
    // can still be used and the component field becomes a real strip ID.
    G4VTouchable* TH = (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();
    G4VSolid* solid = TH->GetSolid();
    G4Trap* Trap = dynamic_cast<G4Trap*>(solid);
    G4Box*  Box  = dynamic_cast<G4Box*>(solid);

    // Readout geometry:
    //   G4Trap -> big uRwell, keep U/V stereo-strip rotation.
    //   G4Box  -> small prototype, use direct local X/Y coordinates.
    bool isXY = (Box != nullptr);

    if (Trap == nullptr && Box == nullptr) {
        cout << "Warning: uRwell sensitive volume is neither G4Trap nor G4Box. Returning original id." << endl;

        if (id.size() > 0) {
            id[id.size() - 1].id_sharing = 1.0;
        }

        return id;
    }

    bool isProto = false;
    if (solid->GetName().find("proto") != std::string::npos) {
        isProto = true;
    }

    if (Trap != nullptr) {
        // Get trapezoid dimensions from the actual G4Trap volume.
        uRwellC.Xhalf_base      = Trap->GetXHalfLength1();
        uRwellC.Xhalf_Largebase = Trap->GetXHalfLength2();
        uRwellC.Yhalf           = Trap->GetYHalfLength1();
        uRwellC.Zhalf           = Trap->GetZHalfLength();
    }
    else {
        // G4Box case: rectangular active area.
        // Geant4 lengths here are in the same CLHEP units used by lxyz.
        uRwellC.Xhalf_base      = Box->GetXHalfLength();
        uRwellC.Xhalf_Largebase = Box->GetXHalfLength();
        uRwellC.Yhalf           = Box->GetYHalfLength();
        uRwellC.Zhalf           = Box->GetZHalfLength();
        /*cout << "DEBUG uRwell G4Box branch used. "
     << "Box Xhalf=" << Box->GetXHalfLength()/mm
     << " Yhalf=" << Box->GetYHalfLength()/mm
     << endl;*/
    }

    double depe = aStep->GetTotalEnergyDeposit();
    double time = aStep->GetPostStepPoint()->GetGlobalTime();

    // If no deposited energy, avoid division by zero later.
    if (depe <= 0) {
        if (id.size() > 0) {
            id[id.size() - 1].id_sharing = 1.0;
        }
        return id;
    }

    // ------------------------------------------------------------
    // STRIP U / X
    // G4Trap: U stereo readout; G4Box: direct local-X readout
    // layer = 1, component = strip number
    // ------------------------------------------------------------
    uRwellC.get_strip_info("strip_u", isProto);
    vector<uRwell_strip_found> multi_hit_u =
        URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto, isXY, "x");
    /*cout << "DEBUG sector=" << id[1].id
     << " layer U, x=" << lxyz.x()/mm
     << " y=" << lxyz.y()/mm
     << " strips found=" << multi_hit_u.size()
     << endl;*/

    //cout << "uRwell U strips found: " << multi_hit_u.size() << endl;

    for (unsigned int h = 0; h < multi_hit_u.size(); h++) {
    /*cout << "  U strip numberID=" << multi_hit_u.at(h).numberID
         << " weight=" << multi_hit_u.at(h).weight
         << endl;*/
        if (multi_hit_u.at(h).numberID<1) continue;

    	/*cout << "  U strip " << h
         << " numberID = " << multi_hit_u.at(h).numberID
         << " weight = " << multi_hit_u.at(h).weight
         << " time = " << multi_hit_u.at(h).time
         << endl;*/
        for (unsigned int j = 0; j < id.size(); j++) {

            identifier this_id;
            this_id.name = id[j].name;
            this_id.rule = id[j].rule;
            this_id.id   = id[j].id;
            this_id.time = id[j].time;

            // U readout plane
            if (this_id.name == "layer") {
                this_id.id = 1;
            }

            // Strip number
            if (this_id.name == "component") {
                this_id.id = multi_hit_u.at(h).numberID;
                this_id.time = multi_hit_u.at(h).time;
            }

            this_id.TimeWindow = id[j].TimeWindow;
            this_id.TrackId    = id[j].TrackId;

            this_id.id_sharing =
                multi_hit_u.at(h).weight / (uRwellC.gain * 1e6 * depe / uRwellC.w_i);

            yid.push_back(this_id);
        }
    }

    // ------------------------------------------------------------
    // STRIP V / Y
    // G4Trap: V stereo readout; G4Box: direct local-Y readout
    // layer = 2, component = strip number
    // ------------------------------------------------------------
    uRwellC.get_strip_info("strip_v", isProto);
    vector<uRwell_strip_found> multi_hit_v =
        URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto, isXY, "y");
    /*cout << "DEBUG uRwell returned IDs: ";
    for (unsigned int ii = 0; ii < id.size(); ii++) {
        cout<<" [id="<<id[ii].id<<", sharing="<<id[ii].id_sharing<<"]";
    }
    cout << endl;*/
    //cout << "uRwell V strips found: " << multi_hit_v.size() << endl;

/*cout << "DEBUG sector=" << id[1].id
     << " layer V, x=" << lxyz.x()/mm
     << " y=" << lxyz.y()/mm
     << " strips found=" << multi_hit_v.size()
     << endl;*/
    for (unsigned int h = 0; h < multi_hit_v.size(); h++) {
  /*  cout << "  V strip numberID=" << multi_hit_v.at(h).numberID
         << " weight=" << multi_hit_v.at(h).weight
         << endl;*/
        if (multi_hit_v.at(h).numberID<1) continue;
	
   /* cout << "  V strip " << h
         << " numberID = " << multi_hit_v.at(h).numberID
         << " weight = " << multi_hit_v.at(h).weight
         << " time = " << multi_hit_v.at(h).time
         << endl;*/

        for (unsigned int j = 0; j < id.size(); j++) {

            identifier this_id;
            this_id.name = id[j].name;
            this_id.rule = id[j].rule;
            this_id.id   = id[j].id;
            this_id.time = id[j].time;

            // V readout plane
            if (this_id.name == "layer") {
                this_id.id = 2;
            }

            // Strip number
            if (this_id.name == "component") {
                this_id.id = multi_hit_v.at(h).numberID;
                this_id.time = multi_hit_v.at(h).time;
            }

            this_id.TimeWindow = id[j].TimeWindow;
            this_id.TrackId    = id[j].TrackId;

            this_id.id_sharing =
                multi_hit_v.at(h).weight / (uRwellC.gain * 1e6 * depe / uRwellC.w_i);

            yid.push_back(this_id);
        }
    }

    return yid;
}


/*vector<identifier> uRwell_HitProcess::processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
    vector<identifier> yid;

    // Local position inside the sensitive cathode gas volume
    G4ThreeVector globalPos = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localPos =
        aStep->GetPreStepPoint()
             ->GetTouchableHandle()
             ->GetHistory()
             ->GetTopTransform()
             .TransformPoint(globalPos);

    double x_mm = localPos.x() / mm;
    double y_mm = localPos.y() / mm;

    // Current Box geometry:
    // dimensions = 10*cm 10*cm zhalf
    // so half-size is 10 cm = 100 mm.
    double half_x_mm = 100.0;
    double half_y_mm = 100.0;

    // LMU pitch is 1 mm.
    double pitch_x_mm = 1.0;
    double pitch_y_mm = 1.0;

    int nstrip_x = int(2.0 * half_x_mm / pitch_x_mm);
    int nstrip_y = int(2.0 * half_y_mm / pitch_y_mm);

    int xstrip = int((x_mm + half_x_mm) / pitch_x_mm) + 1;
    int ystrip = int((y_mm + half_y_mm) / pitch_y_mm) + 1;

    if (xstrip < 1) xstrip = 1;
    if (xstrip > nstrip_x) xstrip = nstrip_x;

    if (ystrip < 1) ystrip = 1;
    if (ystrip > nstrip_y) ystrip = nstrip_y;

    // ------------------------------------------------------------
    // Output logical readout hit 1: X strip
    // layer = 1, component = xstrip
    // ------------------------------------------------------------
    for (unsigned int j = 0; j < id.size(); j++) {
        identifier this_id = id[j];

        if (this_id.name == "layer") {
            this_id.id = 1;
        }

        if (this_id.name == "component") {
            this_id.id = xstrip;
            this_id.id_sharing = 1.0;
            this_id.time = aStep->GetPostStepPoint()->GetGlobalTime();
        }

        yid.push_back(this_id);
    }
    // ------------------------------------------------------------
    // Output logical readout hit 2: Y strip
    // layer = 2, component = ystrip
    // ------------------------------------------------------------
    for (unsigned int j = 0; j < id.size(); j++) {
        identifier this_id = id[j];

        if (this_id.name == "layer") {
            this_id.id = 2;
        }

        if (this_id.name == "component") {
            this_id.id = ystrip;
            this_id.id_sharing = 1.0;
            this_id.time = aStep->GetPostStepPoint()->GetGlobalTime();
        }

        yid.push_back(this_id);
    }

    return yid;
}*/


/*vector<identifier> uRwell_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	
	//uRwellConstants uRwellC;
	
	vector<identifier> yid;
	uRwell_strip URwell_strip;
	// double Lorentz_angle=0;
	G4ThreeVector   xyz    = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector  lxyz    = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); ///< Local Coordinates of interaction
	
	
	
	
	G4VTouchable* TH = (G4VTouchable*) aStep->GetPreStepPoint()->GetTouchable();
	G4Trap *Trap = dynamic_cast<G4Trap*>(TH->GetSolid());
	
	//cout << Trap->GetName()<<endl;

	 bool isProto = false;
	 if(Trap->GetName().find("proto")!=std::string::npos) isProto = true;

	uRwellC.Xhalf_base = Trap->GetXHalfLength1();
	uRwellC.Xhalf_Largebase = Trap->GetXHalfLength2();
	uRwellC.Yhalf = Trap->GetYHalfLength1();
	uRwellC.Zhalf = Trap->GetZHalfLength();
	
	//int sector = id[0].id;
	// int chamber = id[1].id;
	
	double depe = aStep->GetTotalEnergyDeposit();
	double time = aStep->GetPostStepPoint()->GetGlobalTime();
	
	// STRIP U 
    uRwellC.get_strip_info("strip_u", isProto);
    vector<uRwell_strip_found> multi_hit_u = URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto);
    int n_multi_hits_u = multi_hit_u.size();


    for(int h=0; h<n_multi_hits_u; h++){

    		for(int j=0; j<5; j++)
    		{
    			// j=0 region ; j=1 sector; j2 chamber; j3 layer; j4 component

    			identifier this_id;
    			this_id.name       = id[j].name;
    			this_id.rule       = id[j].rule;

    			if(j==0) this_id.id = id[j].id;
    			if(j==1) this_id.id = id[j].id;
    			if(j==2) this_id.id = id[j].id;
    			if(j==3) {
    				this_id.id = 2*id[0].id-1;
    			}
    			this_id.time       = id[j].time;

    			if(j==4){    //J==4 strip ID
    				if(id[2].id>0) {
					    if(multi_hit_u.at(h).numberID ==-15000){
					    	this_id.id  = multi_hit_u.at(h).numberID ;
					    } else {
					    	this_id.id  = multi_hit_u.at(h).numberID + std::accumulate(uRwellC.number_strip_chamber,uRwellC.number_strip_chamber +id[2].id-1,0);
					    }
    				}else this_id.id  = multi_hit_u.at(h).numberID;
    					this_id.time       = multi_hit_u.at(h).time;

    					}

    			this_id.TimeWindow = id[j].TimeWindow;
    			this_id.TrackId    = id[j].TrackId;
    			this_id.id_sharing = multi_hit_u.at(h).weight/(uRwellC.gain*1e6*depe/uRwellC.w_i);

    			yid.push_back(this_id);
    		}

    }

	
    // STRIP V 
    uRwellC.get_strip_info("strip_v", isProto);
	vector<uRwell_strip_found> multi_hit_v = URwell_strip.FindStrip(lxyz, depe, uRwellC, time, isProto);
	int n_multi_hits_v = multi_hit_v.size();

	 for(int h=0; h<n_multi_hits_v; h++){
	 		for(int j=0; j<5; j++)
	 		{
	 			identifier this_id;
				this_id.name       = id[j].name;
				this_id.rule       = id[j].rule;
				if(j==0) this_id.id = id[j].id;
				if(j==1) this_id.id = id[j].id;
				if(j==2) this_id.id = id[j].id;
				if(j==3) this_id.id = 2*id[0].id;
				this_id.time       = id[j].time;

				if(j==4){    //J==4 strip ID
				   if(id[2].id>0) {
					    if(multi_hit_v.at(h).numberID ==-15000){
					    	this_id.id  = multi_hit_v.at(h).numberID ;
					    } else {
					    	this_id.id  = multi_hit_v.at(h).numberID + std::accumulate(uRwellC.number_strip_chamber,uRwellC.number_strip_chamber +id[2].id-1,0);
					    }

				   }else this_id.id  = multi_hit_v.at(h).numberID;
				   	   this_id.time       = multi_hit_v.at(h).time;
				}
				this_id.TimeWindow = id[j].TimeWindow;
				this_id.TrackId    = id[j].TrackId;
				this_id.id_sharing = multi_hit_v.at(h).weight/(uRwellC.gain*1e6*depe/uRwellC.w_i);
				yid.push_back(this_id);
	 	}
	 		}
	
	
	return yid;
	
}*/




// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> uRwell_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;
	
	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);
	
	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)
	
	return noiseHits;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > uRwell_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;
	
	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double uRwell_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}

map< string, vector <int> >  uRwell_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

void uRwell_HitProcess::initWithRunNumber(int runno)
{
	string digiVariation    = gemcOpt.optMap["DIGITIZATION_VARIATION"].args;
	string digiSnapshotTime = gemcOpt.optMap["DIGITIZATION_TIMESTAMP"].args;
	
	if(uRwellC.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		uRwellC = initializeuRwellConstants(runno, digiVariation, digiSnapshotTime, accountForHardwareStatus);
		uRwellC.runNo = runno;
	}
}



// this static function will be loaded first thing by the executable
uRwellConstants uRwell_HitProcess::uRwellC = initializeuRwellConstants(1);











