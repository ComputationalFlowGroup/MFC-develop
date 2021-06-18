#!/bin/bash

# r -> Number of resource sets (RS) per host (node)
# c -> Number of physical cores (CPUs) per RS
# a -> Number of MPI ranks per RS
# g -> Number of GPUs per RS

# -r6 -a7 -c7 -g1 
## Use whole node
## One resource set per GPU (g=1)
## 42/6 = 7 where 42 = number of total cores (c)
## (a == c) usually

# export PGI_ACC_NOTIFY=2
# export NV_ACC_NOTIFY=3
# export NV_ACC_TIME=1
# export NV_ACC_DEBUG=1

# jsrun -r1 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a7 -c7 -g1 ../../src/simulation_code/simulation

# jsrun -r6 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r6 -a7 -c7 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a2 -c2 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a2 -c2 -g1 ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 nsys profile -o output-sys.prof --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

# mpirun -n 1 ../../src/pre_process_code/pre_process
# mpirun -n 1 ../../src/simulation_code/simulation

# mpirun -n 4 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation
# mpirun -n 4 nsys profile -o output-sys.prof --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation


# mpirun -n 1 ncu --launch-skip 2 --launch-count 15 --set full -o output-cu.prof -f   ../../src/simulation_code/simulation
# mpirun -n 1 ncu --launch-skip 2 --launch-count 15 --set full -o output-cu.prof -f --nvtx  ../../src/simulation_code/simulation


## Profile 
# jsrun -r1 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
# jsrun -r1 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

jsrun -r6 -a7 -c7 -g1 ../../src/pre_process_code/pre_process
jsrun -r6 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc,nvtx  --force-overwrite true ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
# jsrun -r1 -a1 -c1 -g1 nvprof --analysis-metrics -o output.nvvp  -f ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ncu --kernel-regex s_weno_alt_177_gpu --launch-skip 3 --launch-count 9 --set full -o output.prof -f  ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ncu --launch-skip 2 --launch-count 5 --set full -o output.prof -f  ../../src/simulation_code/simulation


# jsrun -r1 -a7 -c7 -g1 nsys profile -o output-sys.prof_%q{OMPI_COMM_WORLD_RANK} --stats=true -t openacc  --force-overwrite true ../../src/simulation_code/simulation

# jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
