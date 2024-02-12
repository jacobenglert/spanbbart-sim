#!/bin/bash

array_job1=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 1)
array_job2=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 2)
array_job3=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 3)
array_job4=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 4)
array_job5=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 5)
array_job6=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 6)
array_job7=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 7)
array_job8=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 8)
array_job9=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 9)
array_job10=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 10)
array_job11=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 11)
array_job12=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 12)
array_job13=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 13)
array_job14=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 14)
array_job15=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 15)
array_job16=$(sbatch --parsable HPC/spanbbart-sim-new.slurm 16)

sbatch --depend=afterany:\
$array_job1:\
$array_job2:\
$array_job3:\
$array_job4:\
$array_job5:\
$array_job6:\
$array_job7:\
$array_job8:\
$array_job9:\
$array_job10:\
$array_job11:\
$array_job12:\
$array_job13:\
$array_job14:\
$array_job15:\
$array_job16 \
HPC/spanbbart-sim-combine-new.slurm
