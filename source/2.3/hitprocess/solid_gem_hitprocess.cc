// G4 headers
#include "G4UnitsTable.hh"
#include "G4Poisson.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// gemc headers
#include "solid_gem_hitprocess.h"

map<string, double> solid_gem_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	trueInfos tInfos(aHit);
//	predefined variable Etot, x, y, z, lx, ly, lz, time
	
	int nsteps = aHit->GetPos().size();

	// Get Total Energy deposited
// 	double Etot = 0;
// 	vector<G4double> Edep = aHit->GetEdep();
// 	for(int s=0; s<nsteps; s++) Etot = Etot + Edep[s];
	
	// average global positions of the hit
	
	// Want the first and last local positions
	
// 	double x, y, z;
	double lx_in, ly_in, lz_in;
	double lx_out, ly_out, lz_out;

// 	x = y = z = lx_in = ly_in = lz_in = 0;
	lx_in = ly_in = lz_in = 0;	
	lx_out = ly_out = lz_out = 0;
	vector<G4ThreeVector> pos  = aHit->GetPos();
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	G4ThreeVector p = aHit->GetMom();
	
// 	for(int s=0; s<nsteps; s++) {
// 	    x  = x  +  pos[s].x()/((double) nsteps);
// 	    y  = y  +  pos[s].y()/((double) nsteps);
// 	    z  = z  +  pos[s].z()/((double) nsteps);
// 	}

	/*
	lx_in = Lpos[0].x();
	ly_in = Lpos[0].y();
	lz_in = Lpos[0].z();

	lx_out = Lpos[nsteps-1].x();
	ly_out = Lpos[nsteps-1].y();
	lz_out = Lpos[nsteps-1].z();
	*/

	//  I guess these should be in lab coords
	//  for the digitization software
	lx_in = pos[0].x();
	ly_in = pos[0].y();
	lz_in = pos[0].z();

	lx_out = pos[nsteps-1].x();
	ly_out = pos[nsteps-1].y();
	lz_out = pos[nsteps-1].z();

	// entrance and exit time
	double time_in = 0;
	double time_out = 0;
// 	vector<G4double> times = aHit->GetTime();
	vector<G4double> times = aHit->GetTime();

	time_in = times[0];
	time_out = times[nsteps-1];

	// Energy of the track
	double Ene = aHit->GetE();
	
	dgtz["ETot"] =	tInfos.eTot/MeV;
	dgtz["x"] = tInfos.x/mm;
	dgtz["y"] = tInfos.y/mm;
	dgtz["z"] = tInfos.z/mm;	
	dgtz["lxin"] = lx_in/mm;
	dgtz["lyin"] = ly_in/mm;
	dgtz["lzin"] = lz_in/mm;	
	dgtz["tin"] =  time_in/ns;
	dgtz["lxout"] = lx_out/mm;
	dgtz["lyout"] = ly_out/mm;
	dgtz["lzout"] = lz_out/mm;
	dgtz["tout"] =  time_out/ns;
	dgtz["pid"] =   aHit->GetPID();	
	dgtz["vx"] = aHit->GetVert().getX()/mm;
	dgtz["vy"] = aHit->GetVert().getY()/mm;
	dgtz["vz"] = aHit->GetVert().getZ()/mm;
	dgtz["trE"] = Ene/MeV;
	dgtz["trid"] = aHit->GetTId();
	dgtz["weight"] = 0;	//not used
	dgtz["px"] = p.x()/MeV;
	dgtz["py"] = p.y()/MeV;
	dgtz["pz"] = p.z()/MeV;
	
	dgtz["mpid"]    = (double) aHit->GetmPID();
	dgtz["mtid"]    = (double) aHit->GetmTrackId();
	dgtz["otid"]    = (double) aHit->GetoTrackId();
	dgtz["mvx"]     = aHit->GetmVert().getX();
	dgtz["mvy"]     = aHit->GetmVert().getY();
	dgtz["mvz"]     = aHit->GetmVert().getZ();		
	
	dgtz["id"]  =  identity[0].id;	
	dgtz["hitn"] = hitn;
	
	return dgtz;
}

vector<identifier>  solid_gem_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}


map< string, vector <int> >  solid_gem_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}
