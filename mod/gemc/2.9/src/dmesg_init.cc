/// \file dmesg_init.cc
/// Initialize detector messenger commands.\n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

// G4 headers
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"

// gemc headers
#include "gemcOptions.h"
#include "dmesg_init.h"
#include "string_utilities.h"

/// \fn init_dmesg(goptions)
/// <b> init_dmesg</b>\n\n
/// Initialize detector messenger commands.
vector<string> init_dmesg(goptions gemcOpt)
{
	vector<string> commands;
	
	// Geant4 Process verbosity
	double PHI_VERB        = gemcOpt.optMap["PHY_VERBOSITY"].arg ;
	int   OVERL            = (int) gemcOpt.optMap["CHECK_OVERLAPS"].arg ;
	int   DAWN_N           = (int) gemcOpt.optMap["DAWN_N"].arg ;
	int   SAVE_ALL_MOTHERS = (int) gemcOpt.optMap["SAVE_ALL_MOTHERS"].arg ;
	if (gemcOpt.optMap["SAVE_ALL_ANCESTORS"].arg && (SAVE_ALL_MOTHERS == 0))
	  SAVE_ALL_MOTHERS = 1;
	
	char phi_verb[2];
	sprintf(phi_verb, "%d", (int) PHI_VERB);
	string proc_verb1 = "/process/setVerbose ";
	string proc_verb2 = "/process/verbose ";
	proc_verb1.append(phi_verb);
	proc_verb1.append(" all");
	proc_verb2.append(phi_verb);
	
	
	// Geant4 Tracking Verbosity
	double TRACK_VERB =  gemcOpt.optMap["G4TRACK_VERBOSITY"].arg ;
	char track_verb[2];
	sprintf(track_verb, "%d", (int) TRACK_VERB);
	string tra_verb = "/tracking/verbose ";
	tra_verb.append(track_verb);
	commands.push_back(tra_verb);
	
	// keep trajectories in the container if mother infos are required
	if(SAVE_ALL_MOTHERS)
		commands.push_back("/tracking/storeTrajectory 2");
	else
		commands.push_back("/tracking/storeTrajectory 0");
	
	
	// sets all verbosity to zero
	commands.push_back("/control/verbose 0");
	commands.push_back("/geometry/navigator/verbose 0");
	commands.push_back("/event/verbose 0");
	commands.push_back("/run/verbose 0");
	commands.push_back("/run/particle/verbose 0");
	commands.push_back("/material/verbose 0");
	commands.push_back(proc_verb1.c_str());
	commands.push_back(proc_verb2.c_str());
	commands.push_back("/process/eLoss/verbose 0");
	commands.push_back("/vis/verbose 0");
 	commands.push_back("/vis/viewer/flush");
	commands.push_back("/process/setVerbose 0 all");

	
    // Tell GEANT4 to check for overlaps. This method has changed in GEANT 4.10.0 to check overlap of each
    // volume with its neighbors for about 10000 points on the surface. This is a fairly time consuming procedure.
    // if CHECK_OVERLAPS is set to a number larger than 1, use that many points instead.
	if(OVERL>=1) {
        if (OVERL > 1) {
            commands.push_back("/geometry/test/resolution"+ std::to_string(OVERL));
        }
        commands.push_back("/geometry/test/run");
    }

	if(DAWN_N>0) {
		string DN = "/run/beamOn " + to_string(DAWN_N);

		commands.push_back("/vis/open DAWNFILE");
		commands.push_back("/vis/drawVolume");
		commands.push_back("/vis/scene/add/trajectories rich smooth");
		commands.push_back("/vis/scene/add/hits");
		commands.push_back("/vis/scene/endOfEventAction accumulate -1");
		commands.push_back(DN);
	}
	
	
	
	return commands;
}

/// \fn init_dvmesg(goptions, G4VisManager *VM)
/// <b> init_dvmesg</b>\n\n
/// Initialize detector messenger commands for visualization.
vector<string> init_dvmesg(goptions gemcOpt, G4VisManager *VM)
{
	vector<string> commands;
	commands.push_back("/vis/scene/create gemcLab");
	double use_gui   = gemcOpt.optMap["USE_GUI"].arg;
	int   HIGH_RES  = (int) gemcOpt.optMap["HIGH_RES"].arg ;
	vector<string> wpos =  get_info(gemcOpt.optMap["GUIPOS"].args);
	

	double newxpos = 700 + get_number(wpos[0]);
	string newpos = "+" + stringify(newxpos) + "+" + wpos[1];

	// OpenGL immediate mode
	if(use_gui==1)
		commands.push_back("/vis/open OGLSQt " + gemcOpt.optMap["geometry"].args + newpos);
	// openGL stored mode
	else if(use_gui==2)
		commands.push_back("/vis/open OGLIQt " + gemcOpt.optMap["geometry"].args + newpos);
	

	commands.push_back("/vis/viewer/set/culling coveredDaughters true");


	if(HIGH_RES==1)
	{
		commands.push_back("/vis/viewer/set/background 1 1 1 1");   // 205, 230, 251 = CD, E6, FA
		commands.push_back("/vis/viewer/set/autoRefresh 1");
		commands.push_back("/vis/scene/add/trajectories rich smooth");
		commands.push_back("/vis/scene/add/hits");
		// This rendering can be problematic for some operation solid.
		// Leaving this out for now, and letting G4 do what's more comfy.
		// commands.push_back("/vis/viewer/set/lineSegmentsPerCircle 50");
		commands.push_back("/vis/scene/endOfEventAction accumulate -1");

		// increasing visualisation storing
		commands.push_back("/vis/ogl/set/displayListLimit 500000");
	}
	else
	{
		commands.push_back("/vis/scene/add/trajectories");
		commands.push_back("/vis/scene/add/hits");
		commands.push_back("/vis/scene/endOfEventAction accumulate -1");
	}
	
	

	
	// tracks colors
	// these steps come from the user manual chapter Trajectory Drawing Models
	// Create and configure a drawByParticleID model named  gemcColorIDModel
	// positive particles: RED
	// negative particles: GREEN
	// specific colors below
	G4TrajectoryDrawByParticleID* gemcColorIDModel = new G4TrajectoryDrawByParticleID("gemcColorIDModel");
	gemcColorIDModel->SetDefault("gray");
	gemcColorIDModel->Set("neutron", "black");
	gemcColorIDModel->Set("gamma",   "blue");
	gemcColorIDModel->Set("e-",      "red");
        gemcColorIDModel->Set("e+",      "yellow");        
	gemcColorIDModel->Set("pi+",     "magenta");
	gemcColorIDModel->Set("pi-",     "magenta");
        gemcColorIDModel->Set("mu+",     "brown");
        gemcColorIDModel->Set("mu-",     "brown");
	gemcColorIDModel->Set("proton",  G4Colour(0.95, 0.6, 0.3));  // orange
        gemcColorIDModel->Set("opticalphoton", "white");  
        
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	for(int i=0; i<particleTable->entries(); i++)
	{
		string pname  = particleTable->GetParticleName(i);
		double charge = particleTable->FindParticle(pname)->GetPDGCharge();
		string pcolor ;
		if(pname !=  "neutron" &&
		   pname !=  "gamma"   &&
		   pname !=  "e-"      &&
                   pname !=  "e+"      &&		   
		   pname !=  "pi+"     &&
		   pname !=  "pi-"     &&
                   pname !=  "mu+"     &&
                   pname !=  "mu-"     &&
		   pname !=  "proton"  &&
                   pname !=  "opticalphoton"     )
		{
			if(charge>0)
				gemcColorIDModel->Set(pname,  "cyan");
			if(charge==0)
				gemcColorIDModel->Set(pname,  "gray");
			if(charge<0)
				gemcColorIDModel->Set(pname,  "green");
		}
	}
	
	VM->RegisterModel(gemcColorIDModel);
	VM->SelectTrajectoryModel(gemcColorIDModel->Name());
	
	// Draw by Process - interesting! But not working yet?
	// 	commands.push_back("/vis/modeling/trajectories/create/drawByAttribute gemcColorAttributeModel");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/setAttribute CPN");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/verbose true");
	//
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/addValue brem_key       eBrem");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/addValue annihil_key annihil");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/addValue decay_key Decay");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/addValue muIon_key muIoni");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/addValue eIon_key  eIoni");
	//
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/brem_key/setLineColour    yellow");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/annihil_key/setLineColour    yellow");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/decay_key/setLineColour    yellow");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/eIon_key/setLineColour    yellow");
	// 	commands.push_back("/vis/modeling/trajectories/gemcColorAttributeModel/muIon_key/setLineColour    yellow");
	
	return commands;
}






















