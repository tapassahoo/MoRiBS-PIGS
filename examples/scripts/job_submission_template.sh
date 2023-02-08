#!/bin/bash
#SBATCH --job-name=itcfn2p15000
#SBATCH --output=itcfn2p15000.log
#SBATCH --error=itcfn2p15000.err
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

var_job_type=analysis
var_molecule=2
var_parameter=0.2
var_preskip=15000

echo "job_type       = $var_job_type"
echo "numb_molecule  = $var_molecule"
echo "numb_parameter = $var_parameter"
echo "numb_preskip   = $var_preskip"

cd /home/tapas/MoRiBS-PIGS/examples/scripts/
cp generic_execution_moribs_driver_linear_molecule_itcf.py temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py
sed -i "s/\<JOB_TYPE\>/$var_job_type/" temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py
sed -i "s/\<NUMB_MOLECULE\>/$var_molecule/" temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py
sed -i "s/\<PARAMETER_VALUE\>/$var_parameter/" temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py
sed -i "s/\<NUMB_PRESKIP\>/$var_preskip/" temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py
python temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py

date
echo -e "\n"
echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
