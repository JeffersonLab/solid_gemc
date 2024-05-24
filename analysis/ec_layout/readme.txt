= Intro and howto ========================
From the drawing and tables by engineers, we got the EC hexagon module center coordinate in lab frame and save as txt file

EC_map.C script tool
1. root 'EC_map.C("find",100,100)'
   get the fired module and its surrounded 6 module ID with hit position hit_x hit_y in cm
2. root EC_map.C("convert") 
   convert the EC coordinate map "map_FAEC_ANL_20130628.txt" into GEMC parameter file "map_FAEC_ANL_20130628_PVDIS_parameters.txt"
   
ec_layout.h lib tool can be used as follows and an example in analysis/ml/analysis_aiml.C
*********************
#include "../ec_layout/ec_layout.h"
LoadEC_map("../ec_layout");
int ecal_hitblockID=GetECALBlock_id(hit_x,hit_y); //get ec id from hit position
const int Nneighbor=19; //hexagon has 6+1 or 18+1 layout
Int_t cluster_edep_blockid[Nneighbor]={0};      
int block_count = GetECALCluserAll(Nneighbor,ecal_hitblockID,cluster_edep_blockid);             
*********************

= details ===============================

the SoLID lab frame
*********************
The solenoid coil central plane and solenoid central line crosses at the origin.
electron beam and solenoid central line are along +z axis
y axis is vertical and +y pointing up relative to the ground
x axis is horizontal and +x pointing left relative to the electron beam
*********************
refer to https://solid.jlab.org/wiki/index.php/Solid_Software#Coordinate_System

Imaging we are riding along the beam and facing EC from front
EC module ID consist of idx,idy,id,sector
idx and idy are given by column and row number of the module
Their values increase from bottom to top for idy, from right to left for idx
(idx=1,idy=1) module is the most bottom right module and its coordinate x<0 and y<0
In the script, module center coordinate in lab frame is recorded in x[idy][idx], y[idy][idx]
id is the total id of the module, its value increase by following idx first, then idy
sector is the sector number, it matches baffle for FAEC

The module layout is for PVDIS FAEC
the range are
sector: 1-30 
id:1-1748
idy: 1-54 rows
idx direction: 1-46 coll. in maximum row
status can be 1 (default) or 0, which can be used to enable or disable a particular module

SIDIS FAEC and LAEC are using the PVDIS FAEC module layout with certain modules removed when not in the right radius range.
To remove module, we can just use radius from x y info instead of keyword status

One hidden information in ANL layout is the angle of hexagon module.
It has two side parallel to y-axis, no side parallel to x-axis

= log ====================================================

written by Yuxiang Zhao (yxzhao@jlab.org)  2014/08

Updated by Zhiwen Zhao  (zwzhao@jlab.org) and Ye Tian (tianye@jlab.org)
