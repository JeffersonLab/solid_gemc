#!/bin/csh -f

source /jlab/2.5/ce/jlab.csh
setenv GEMC /group/solid/apps/gemc/2.9
setenv SoLID_GEMC /solid_gemc
setenv LD_LIBRARY_PATH ${GEMC}:${LD_LIBRARY_PATH}
setenv PATH ${SoLID_GEMC}/source/${GEMC_VERSION}:${PATH}

#run simulation
# solid_gemc solid.gcard -BEAM_P="e-,5.5*GeV,20*deg,0*deg" -SPREAD_P="5.5*GeV,15*deg,180*deg" -BEAM_V="(0,0,-350)cm" -SPREAD_V="(0,25)cm" -N=1e4 | cat >> out.log

solid_gemc solid.gcard -INPUT_GEN_FILE="LUND,gen.lund" -N=1e4 | cat >> out.log

# solid_gemc solid.gcard -FILTER_HITS=1 -BEAM_P="e-,11*GeV,0*deg,0*deg" -SPREAD_P="0*GeV,0*deg,0*deg" -BEAM_V="(0,0,-400)cm" -SPREAD_V="(0.25,0)cm"  -N=5e7 | cat >> out.log

#convert evio to root
#need more option for banks other than flux
# which  evio2root
evio2root -INPUTF=out.evio -B="/solid_gemc/geometry/ec_segmented_moved/solid_SIDIS_ec_forwardangle /solid_gemc/geometry/ec_segmented_moved/solid_SIDIS_ec_largeangle /solid_gemc/geometry/hgc_moved/solid_SIDIS_hgc /solid_gemc/geometry/lgc_moved/lg_cherenkov  /solid_gemc/geometry/gem_moved/solid_SIDIS_gem /solid_gemc/geometry/mrpc_moved/solid_SIDIS_mrpc_forwardangle /solid_gemc/geometry/spd_moved/solid_SIDIS_spd_forwardangle /solid_gemc/geometry/spd_moved/solid_SIDIS_spd_largeangle" -R="flux" | cat >> out.log

# cp /solid_gemc/analysis/SoLIDFileReduce.C ./
# root -b -q 'SoLIDFileReduce.C+("SIDIS","out.root","2")'

# mv -f out.root_skim out.root


