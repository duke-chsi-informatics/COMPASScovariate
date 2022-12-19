#!/bin/bash
#
#SBATCH -p chsi
#SBATCH -A chsi
#SBATCH --time=08:00:00
#SBATCH --mem=1G
#SBATCH -c 4

export SINGULARITY_IMAGE="oras://gitlab-registry.oit.duke.edu/chsi-informatics/containers/compasscovariate-singularity:v003"
export BIND_ARGS="--bind /work:/work"

export ITER=1000
export REP=1

singularity exec ${BIND_ARGS} ${SINGULARITY_IMAGE} Rscript -e "rmarkdown::render('$HOME/project_repos/COMPASScovariate/incubator/simu_cov.Rmd')"
