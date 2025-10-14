#!/bin/tcsh
uname -a | cat > out.log
more /etc/redhat-release | cat >> out.log
lscpu | cat >> out.log
more /proc/meminfo | cat >> out.log
pwd | cat >> out.log

setenv location /scigroup/cvmfs/halla/solid/soft
#run generator
#singularity exec -B /group:/group -B /u:/u -B /w/work:/work -B /w:/w -B /cache:/cache -B /volatile:/volatile -B /lustre:/lustre -B $location/solidevgen_tag1:/evgen  $location/container/jeffersonlab_solidevgen_tag1_latest.sif ./evgen.sh commit0acacfe_20230908
# run solid_gemc 
# singularity exec -B ${PWD}:/mywork -B /group:/group -B /u:/u -B /w/work:/work -B /w:/w -B /cache:/cache -B /volatile:/volatile -B /lustre:/lustre -B $location/jlabce_tag2.5/solid_gemc/commit4dbc836_20220531:/solid_gemc -B $location/field:/field $location/container/jeffersonlab_jlabce_tag2.5_digest:sha256:9b9a9ec8c793035d5bfe6651150b54ac298f5ad17dca490a8039c530d0302008_20220413_s3.9.5.sif ./do_it_all.sh
singularity exec -B ${PWD}:/mywork -B /group:/group -B /u:/u -B /w/work:/work -B /w:/w -B /cache:/cache -B /volatile:/volatile -B /lustre:/lustre -B /group/solid/solid_github/JeffersonLab/solid_gemc_test:/solid_gemc -B $location/field:/field $location/container/jeffersonlab_jlabce_tag2.5_digest:sha256:9b9a9ec8c793035d5bfe6651150b54ac298f5ad17dca490a8039c530d0302008_20220413_s3.9.5.sif ./do_it_all.sh

