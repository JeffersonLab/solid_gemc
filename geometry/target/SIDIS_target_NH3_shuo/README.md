# For the SoLID target part
  This SoLID NH3 target part definition. The stl files are imported from the CAD file. The arrangement is in cad.gxml. 
  I converted some parts of the CAD file to geant4 volumes in solid_SIDIS_target_NH3_geometry.pl
. The converted volume stl files are then deleted. 
##Volumes converted
```bash
-make_target_fiels();
-make_scattering_chamber();
-make_scattering_windows();
-make_target();#NH3
-make_target_endcaps();#Al
-make_target_LHe();#liquid He
-make_target_tail_nose();
-make_target_4Kshield();#Al
-make_target_LN2shield();
-make_target_steel();#stainless steel 
-make_target_coil();#Cu coils
-make_target_coil_box();#Al coil holder
-make_target_coil_box_2();#Al coil holder
-make_target_coil_lid();#Al coil holder lid
```
##To use
First build the geometry
```bash
./solid_SIDIS_target_NH3_geometry.pl config_solid_SIDIS_target_NH3.dat
```
Then run the gcard file
```bash
solid_gemc nh3.gcard
```
