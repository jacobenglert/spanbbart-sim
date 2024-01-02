#!/bin/bash

array_job0=$(sbatch --parsable HPC/spanbbart-sim.slurm)

sbatch --depend=afterany:$array_job0 HPC/spanbbart-sim-combine.slurm
