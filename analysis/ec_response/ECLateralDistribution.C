/*
 * ECLateralDistribution.C
 *
 *  Created on: Nov 6, 2013
 *      Author: jinhuang
 */

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>

enum enu_Calo_Config
{
	kPVDIS, kSIDIS_FORWARD, kSIDIS_LARGE
};
enum enu_Particle_Type
{
	kElectron, kPion
};

TH3* h3s[10][10] =
{ NULL };
TH2D* h2s[10][10][100] =
{ NULL };
TFile * f = NULL;

//! Get 2D distribution for EC energy deposition
//! \param[in]	cfg		Calorimeter configuration
//!	\param[in]	part	Particle type
//! \param[in]	TotalEnergyBeforeImpact		Total energy of the particle immediate
//! \output				Energy distribution in the R-phi plane in the unit of mm. (0,0) is the impact location at the surface of calorimeter. z axis is relative energy per R-phi bin, which is normalized to 1.
TH2D *
ECLateralDistribution(const enu_Calo_Config cfg, const enu_Particle_Type part,
        const double TotalEnergyBeforeImpact)
{
	f = new TFile("ECLateralDistribution.root");
	assert(f);

	if ( !h3s[cfg][part] )
	{
		h3s[cfg][part] = f->GetObjectChecked(
		        Form("hHit_XY_VtxE_%d_%d", cfg, part), "TH3");
		assert(h3s[cfg][part]);
	}

	TH3 * h3 = h3s[cfg][part];
	assert(h3);

	const int bin = h3->GetZaxis()->FindBin(TotalEnergyBeforeImpact);

	TH2D * h2 = h2s[cfg][part][bin];
	if ( !h2 )
	{
		h3->GetZaxis()->SetRange(bin, bin);
		h2 = (TH2D *) h3->Project3D(Form("xy_%d", bin));
		assert(h2);
		h2->Scale(1. / h2->GetSum());
		h2->SetTitle("Relative energy deposition;phi*R - EC surface entrance point (mm);R - EC surface entrance point (mm);Relative energy deposition");
		h2s[cfg][part][bin] = h2;
	}
	assert(h2);

	return h2;

}

//! make the data file from simulation
void ExtractHistogram(void)
{
	TFile * f = new TFile("ECLateralDistribution.root", "recreate");
	assert(f);

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.SIDIS_Forward_Run.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone = (TH3*) h->Clone(
	        Form("hHit_XY_VtxE_%d_%d", kSIDIS_FORWARD, kELectron));
	hClone -> Write();

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.SIDIS_Forward_RunPi.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone = (TH3*) h->Clone(
	        Form("hHit_XY_VtxE_%d_%d", kSIDIS_FORWARD, kPion));
	hClone -> Write();

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.SIDIS_Large_Run.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone = (TH3*) h->Clone(
	        Form("hHit_XY_VtxE_%d_%d", kSIDIS_LARGE, kELectron));
	hClone -> Write();

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.SIDIS_Large_RunPi.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone = (TH3*) h->Clone(
	        Form("hHit_XY_VtxE_%d_%d", kSIDIS_LARGE, kPion));
	hClone -> Write();

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.PVDIS_Run.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone =
	        (TH3*) h->Clone(Form("hHit_XY_VtxE_%d_%d", kPVDIS, kELectron));
	hClone -> Write();

	TFile * f1 = new TFile(
	        "script/Lead2X0PbBlock_HalfcmSpt_Hex.1.PVDIS_RunPi.root");
	assert(f1);

	TH3 * h = f1->GetObjectChecked("hHit_XY_VtxE", "TH3");
	assert(h);
	f->cd();
	TH3* hClone = (TH3*) h->Clone(Form("hHit_XY_VtxE_%d_%d", kPVDIS, kPion));
	hClone -> Write();

}
