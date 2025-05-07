/*
 * FastResponse.C
 *
 *  Created on: Aug 15, 2013
 *      Author: jinhuang
 */

#include <iostream>
#include <cassert>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TMath.h>

using namespace std;

//! Get a quick response model for Shashlik calorimeter for a particle before entering shower + preshwoer
//! \param[in] particle names, pi, p, e, gamma
//! \param[in] section PS or S
//! \param[in] depositiontype scint or whole
//! \param[in] E_Kine kinematic in GeV
//! \return Energy deposition in shower. Stocastically generated shower response
double FastResponse(const TString & particle,const TString & section, const TString & depositiontype, const double E_Kine)
{
	// from G4 simulation of Shashlik configuration of 0.5mm Pb+1.5mm Scint.
	//shower calibration by a constant
	static const double E_calibration = 3.88271;
	//preshower calibration by a formula
	//38.9 MeV + 2.549 * ps_raw

	const char * data_file = NULL;

	double mass = 0;

	if (particle == "e")
	{
		data_file =
		        "Elec.10.11.12.root_DrawLowEBgd_EDepShowerOverE_Dist_Read.root";
		mass = 0.510998928e-3;
	}
	else if (particle == "pi")
	{
		data_file =
		        "Pion.run10x.root_DrawLowEBgd_EDepShowerOverE_Dist_Read.root";
		mass = 139.570e-3;
	}
	else if (particle == "p")
	{
		data_file =
		        "Proton.run10x.root_DrawLowEBgd_EDepShowerOverE_Dist_Read.root";
		mass = 938.27204e-3;
	}
	else if (particle == "gamma")
	{
		data_file =
		        "Gamma.run10x.root_DrawLowEBgd_EDepShowerOverE_Dist_Read.root";
		mass = 0;
	}
	else
	{
		cout << "FastResponse - Error - " << "unsupported particle " << particle
		        << endl;

	}

	TFile * file = (TFile *) gROOT->GetListOfFiles()->FindObject(data_file);

	if (!file)
	{
		file = new TFile(data_file);

		assert(file);

		cout <<"FastResponse - INFO - openned "<<data_file<<endl;
	}

	TH2F * h2;
	if (section == "PS") h2 = (TH2F *) file->GetObjectChecked("EDepPreshower_Dist", "TH2F");	
	else if (section == "S") h2 = (TH2F *) file->GetObjectChecked("EDepShowerOverE_Dist", "TH2F");
	else cout << "FastResponse - Error - " << "unsupported section " << section  << endl;
	
	assert(h2);

	const double log_Ek = TMath::Log10(E_Kine);
	const double E_total = TMath::Sqrt(E_Kine * E_Kine + mass * mass);
	const int bin = h2->GetYaxis()->FindBin(log_Ek);

	TH1D* h1 = h2->ProjectionX("proj", bin, bin);

	double scintillator_Edep=0.,whole_Edep=0.;		
	if (section == "PS") {
	  scintillator_Edep = h1->GetRandom()/1e3;  //result in GeV
	  whole_Edep = (38.9+2.549*h1->GetRandom())/1e3;  //result in GeV	  	  
	}
	else if (section == "S") {
	  scintillator_Edep = h1->GetRandom() * E_total; //result in GeV
	  whole_Edep = h1->GetRandom() * E_total * E_calibration; //result in GeV
	}
	else cout << "FastResponse - Error - " << "unsupported section " << section  << endl;

	double Edep=0;
	if (depositiontype == "scint") Edep=scintillator_Edep;
	else if (depositiontype == "whole") Edep=whole_Edep;
	else cout << "FastResponse - Error - " << "unsupported depositiontype " << depositiontype  << endl;	

	return Edep;

}

