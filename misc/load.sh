#!/usr/bin/env bash

echo -e "\nPlease run this file with \"source\": \"source load.sh\"\n".

echo "Which computer would you like to load submodules for?"
echo " - Summit     (s)"
echo " - Bridges2   (b)"
echo " - Ascent     (a)"
echo " - Wombat     (w)"
echo " - Richardson (r)"
echo -n "(s/b/a/w/r): "
read response

if [ "$response" == "s" ]; then # For Summit
    module purge
    module load lsf-tools/2.0
    module load darshan-runtime/3.3.1-lite
    module load DefApps
    module load spectrum-mpi/10.4.0.3-20210112
    module load hsi/5.0.2.p5
    module load xalt/1.2.1
    module load nvhpc/21.9
    module load cuda/11.4.2
    module list
elif [ "$response" == "b" ]; then # Bridges2
    module purge
    module load openmpi/4.0.5-nvhpc21.7
    module load cuda/11.1.1
    module list
elif [ "$response" == "a" ]; then # For Ascent
    module purge
    module load nvhpc/21.9
    module load spectrum-mpi
    module load cuda/11.2.0
    module load nsight-compute
    module load nsight-systems
    module list
elif [ "$response" == "r" ]; then # Richardson
    module purge
    module load gcc/9.3.0
    module load openmpi-2.0/gcc-9.3.0
    module load python/3.7
    module list
elif [ "$response" == "w" ]; then # For Wombat
    module purge
    module load python/3.9.9
    module load /sw/wombat/Nvidia_HPC_SDK/modulefiles/nvhpc/22.1
    module list
else
    echo "Error: Requested system is not supported"
fi
