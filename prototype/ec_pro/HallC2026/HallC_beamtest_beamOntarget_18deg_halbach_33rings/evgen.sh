#!/bin/csh -f

# echo "ok"
source /w/halla-scshelf2102/solid/tianye/solid_prelead2x0/solidevgen_tag1/evgen_inclusive_e/$1/setup
ln -s /w/halla-scshelf2102/solid/tianye/solid_prelead2x0/solidevgen_tag1/evgen_inclusive_e/$1/dat 
/w/halla-scshelf2102/solid/tianye/solid_prelead2x0/solidevgen_tag1/evgen_inclusive_e/$1/evgen_inclusive_e input.dat 1  
