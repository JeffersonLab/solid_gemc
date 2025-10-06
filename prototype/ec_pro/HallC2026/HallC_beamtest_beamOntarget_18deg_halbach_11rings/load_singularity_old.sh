#!/bin/tcsh
uname -a | cat > out.log
more /etc/redhat-release | cat >> out.log
lscpu | cat >> out.log
more /proc/meminfo | cat >> out.log
pwd | cat >> out.log

#may need to slurm job
#source /etc/profile.d/modules.csh  
source /u/home/tianye/cshrc_zwzhao
module load singularity/3.9.5
which singularity | cat >> out.log

setenv location /scigroup/cvmfs/halla/solid/soft
# setenv location /cvmfs/oasis.opensciencegrid.org/jlab/halla/solid/soft

#run generator
#singularity exec -B  $location/solidevgen_tag1:/evgen  $location/container/jeffersonlab_solidevgen_tag1_latest.sif ./evgen.sh commit517d0c6_20220527

# run solid_gemc 
singularity exec -B $location/jlabce_tag2.5/solid_gemc/commit4dbc836_20220531:/solid_gemc -B $location/field:/field $location/container/jeffersonlab_jlabce_tag2.5_digest:sha256:9b9a9ec8c793035d5bfe6651150b54ac298f5ad17dca490a8039c530d0302008_20220413_s3.9.5.sif ./do_it_all.sh

# singularity shell -s /bin/tcsh -B $location/jlabce_tag2.5/solid_gemc/commit4dbc836_20220531:/solid_gemc -B $location/field:/field $location/container/jeffersonlab_jlabce_tag2.5_digest:sha256:9b9a9ec8c793035d5bfe6651150b54ac298f5ad17dca490a8039c530d0302008_20220413_s3.9.5.sif 

