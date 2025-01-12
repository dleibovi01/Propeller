#!bin/bash
export TIME="02:00:00"
export PARTITION="gg"
export NNODES=1


export job_name="nvr_lpr_misc-openfoam.dev"
export gpus_per_node=1
readonly _cont_image="gitlab-master.nvidia.com/modulus/solvers/gpu-openfoam-docker-build/amd64:cuda-12.4-0.1"
readonly _cont_name="openfoam-dev"
readonly _cont_mounts="/lustre/fsw/portfolios/nvr/users/dleibovici/git_repos/gpu-openfoam-docker-build/:/code:rw"

# srun -A coreai_climate_earth2 -p ${PARTITION} -N ${NNODES} -t ${TIME} \
#     --job-name="${job_name}" \
#     --ntasks-per-node=16 \
#     --gpus-per-node=${gpus_per_node} \
#     --container-image="${_container_image}" \
#     --container-name="${_container_name}" \
#     --container-mounts="${_container_mounts}" \
#     --container-workdir /code \
#     --pty /bin/bash -i &

export PYTHONPATH="/code"


srun \ 
  -p ${PARTITION} \
  -N 1 \
  -n 72 \
  --container="${_cont_image}" \
#  --container-mounts="${_cont_mounts}" \
#   --container-workdir=/code \
#  --container-name="${_cont_name}" \
#  --no-container-entrypoint \
#  --mpi=pmix \
  -J "oth24023" 
#  --time=240 \
#  --pty bash

