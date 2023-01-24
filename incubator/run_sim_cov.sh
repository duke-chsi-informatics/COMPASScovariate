#!/bin/bash
#
#SBATCH -p chsi
#SBATCH -A chsi
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -c 16

#Singularity container setup
export SINGULARITY_IMAGE="oras://gitlab-registry.oit.duke.edu/chsi-informatics/containers/compasscovariate-singularity:v003"
export BIND_ARGS="--bind /work:/work"

#Simulation arguments
export ITER=100000
export REP=2
export SIM_REP=5

#Run the simulation
singularity exec ${BIND_ARGS} ${SINGULARITY_IMAGE} Rscript -e "rmarkdown::render('$HOME/project_repos/COMPASScovariate/incubator/simu_cov.Rmd')"

#Copy the knit HTML to /work
cp "${HOME}/project_repos/COMPASScovariate/incubator/simu_cov.html" "/work/${USER}/COMPASScovariate-sim/output"
