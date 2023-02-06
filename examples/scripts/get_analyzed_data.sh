#!/bin/bash
echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo -e "\n"
date
var_job_type=$1
var_molecule=$2
var_parameter=$3
var_preskip=$4
pwd
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
python temp_execution_moribs_driver_linear_molecule_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_analysis.py > ${HOME}/final-pigs-outputs-for-plotting/moribs_${var_job_type}_n${var_molecule}_beta${var_parameter}_nskip${var_preskip}_itcf.log &
date
echo -e "\n"
echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

