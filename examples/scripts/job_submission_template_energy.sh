#!/bin/bash
#SBATCH --job-name=energy_n2_p0
#SBATCH --output=energy_n2_p0.log
#SBATCH --error=energy_n2_p0.err
#SBATCH --time=00-0:30
#SBATCH --account=rrg-pnroy
#SBATCH --constraint=broadwell
#SBATCH --mem-per-cpu=1024
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
echo $SLURM_JOB_ID

echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo -e "\n"
date
set -xv

job_type=analysis
numb_molecule=2
parameter_value=0.1
numb_block=10000
numb_pass=500
preskip_value=0

echo "job_type        = $job_type"
echo "numb_molecule   = $numb_molecule"
echo "parameter_value = $parameter_value"
echo "numb_block      = $numb_block"
echo "numb_pass       = $numb_pass"
echo "numb_preskip    = $preskip_value"

cd /home/tapas/MoRiBS-PIGS/examples/scripts/
cp generic_execution_moribs_driver_linear_molecule.py temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_JOB_TYPE\>/$job_type/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_NUMB_MOLECULE\>/$numb_molecule/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_PARAMETER_VALUE\>/$parameter_value/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_NUMB_BLOCK\>/$numb_block/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_NUMB_PASS\>/$numb_pass/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
sed -i "s/\<INPUT_NUMB_PRESKIP\>/$preskip_value/" temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py
python temp_execution_moribs_driver_linear_molecule_n${numb_molecule}HF_beta${parameter_value}inverse_kelvin_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py

date
echo -e "\n"
echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
