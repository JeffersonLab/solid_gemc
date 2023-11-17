# For the SoLID SIDIS Transverse pol NH3 target using magnet from Scientific Magnet
The magnet was delivered to jlab in late 2021. 
we use the magnet cad file from jlab target group and also refer g2p target to make the target in simulation

build the geometry
```
./solid_SIDIS_target_NH3.pl config_solid_SIDIS_target_NH3.dat

./solid_SIDIS_target_NH3.pl config_solid_SIDIS_target_NH3.dat cad
```
run the gcard file to show graphic
```
solid_gemc nh3.gcard
```

geometry details
* All geometries. 
```
-make_target_field();#a box of air containing the entire target
-make_scattering_chamber();#Tube
-make_scattering_windows();
-make_target();#NH3_He
-make_target_endcaps();#Al
-make_target_cell();#Target cell mounted on the ladder
-make_target_LHe();#liquid He
-make_target_LHe_shield();#liquid He container
-make_target_steel();#stainless steel 
-make_target_coil();#Cu coils
-make_target_coil_box();#Al coil holder
-make_target_coil_box_2();#Al coil holder
-make_target_coil_lid();#Al coil holder lid
-make_40K_shield();#the out tube tank
_make_magnet_support;#defines the 25deg acceptance
```

2023/10/25 first checking, Shuo Jia and Zhiwen Zhao
2023/11/08 updated the target geometry, removed the magnet support cad file. Shuo Jia
