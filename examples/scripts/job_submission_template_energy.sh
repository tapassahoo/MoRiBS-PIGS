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

echo "#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#"
date
echo ""
#set -xv

extra_path="academic-project/"
job_type=submission
simulation_type=PIMC
numb_molecule=2
parameter_value=0.1
numb_block=10
numb_pass=10
preskip_value=0

echo "job_type        = $job_type"
echo "simulation_type = $simulation_type"
echo "numb_molecule   = $numb_molecule"
echo "parameter_value = $parameter_value"
echo "numb_block      = $numb_block"
echo "numb_pass       = $numb_pass"
echo "numb_preskip    = $preskip_value"

new_file="temp_file_for_execution_of_moribs_${simulation_type}_driver_for_linear_molecule_${numb_molecule}HF_beta${parameter_value}kelvin_inverse_mc_blocks${numb_block}_mc_passes${numb_pass}_nskip${preskip_value}_${job_type}.py";

echo ""
echo "Name of the script file for the submission is - $new_file";
echo ""

cd ${HOME}/${extra_path}MoRiBS-PIGS/examples/scripts/
cp generic_execution_moribs_driver_linear_molecule.py  $new_file

sed -i "s/\<INPUT_JOB_TYPE\>/$job_type/"               $new_file
sed -i "s/\<INPUT_SIMULATION_TYPE\>/$simulation_type/" $new_file
sed -i "s/\<INPUT_NUMB_MOLECULE\>/$numb_molecule/"     $new_file
sed -i "s/\<INPUT_PARAMETER_VALUE\>/$parameter_value/" $new_file
sed -i "s/\<INPUT_NUMB_BLOCK\>/$numb_block/"           $new_file
sed -i "s/\<INPUT_NUMB_PASS\>/$numb_pass/"             $new_file
sed -i "s/\<INPUT_NUMB_PRESKIP\>/$preskip_value/"      $new_file

python $new_file

echo ""
date
echo "#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#"
