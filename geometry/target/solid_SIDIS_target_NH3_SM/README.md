# For the SoLID SIDIS Transverse pol NH3 target using magnet from Scientific Magnet
The magnet was delivered to jlab in late 2021. 
we use the magnet cad file from jlab target group and also refer g2p target to make the target in simulation

build the geometry
```
./solid_SIDIS_target_NH3.pl config_solid_SIDIS_target_NH3.dat
```
run the gcard file to show graphic
```
solid_gemc nh3.gcard
```

geometry details
* the magnet coil support is loaded from cad file in cad.gxml
* some are from cad file and some are from old solid_SIDIS_target_NH3_geometry.pl
```
-make_target_field();#a box of air containing the entire target
-make_scattering_chamber();#Tube
-make_scattering_windows();
-make_target();#NH3_He
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

2013/10/25 first checking, Shuo Jia and Zhiwen Zhao