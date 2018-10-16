Here are few instructions that will help a user to modify script files before submitting jobs in queue. The source codes can be downloaded directly from github by typing the following line on the terminal:

```git clone git@github.com:tapassahoo/MoRiBS-PIGS.git```

> It is important to note that PotFunc() function in mc_estim.cc file includes analytic potential like dipole-dipole interaction and the unit of energy is Kelvin.

First, the author must read README file placed in MoRiBS-PIGS/ directory and follow the instructions.

In the source directory, there are many Makefiles. Makefile-PIMC and Makefile-PIGS are the makefiles that a user needs to compile the source codes for finite temperature (PIMC) and ground state (PIGS) canculations, respectively. To compile the source codes, first copy Makefile-PIMC of Makefile-PIGS to Makefile and use the following command:

```
make clean
make
```

But the user does not need to compile the source codes manually if the user wish to submit jobs by script files. The script files are in dir: MoRiBS-PIGS/examples/scripts). There are three python scripts:

- [x] script_submission_analysis_MoRiBS.py

- [x] support.py

- [x] inputFile.py

The user are suggested to make the following modifications in the scripts before running MoRiBs successfully:

1. In **script_submission_analysis_MoRiBS.py**

   - If user wish to run MoRiBs in graham.computecanada.ca, just replace "NameOfServer = "nlogn"" by "NameOfServer = "graham"". "NameOfServer = "nlogn"" when jobs will be submitted in feynman or nlogn server.

   - If the user wish to include cage potential, he/she should use "status_cagepot = True", otherwise, "status_cagepot = False".

   - Keep the same directotory-tree as as the developer used - /home/user_name/source_dir/input_dir. The user may change the names of the directories. As for example, the developer used

```
user_name           = "tapas"
source_dir          = "Moribs-pigs/MoRiBS-PIMC/"                    #Path of the source directory#
out_dir             = "nonlinear-molecule/"                         #This directory will automatically be created in /work or in /scratch if it does not exits.      
input_dir           = "examples/nonlinear-molecule/"                #Where all the input and scripts are
final_results_path  = "/home/"+user_name+"/ResultsOf"+TypeCal+"/"   #Where all the final results will be stored after analyzing the MoRiBs outputs.
```

   - but the user may change these as

```
user_name           = "user_name" excluding "/"
source_dir          = "MoRiBS-PIMC/"
out_dir             = "PIMC-H2O/"
input_dir           = "INPUT/"
final_results_path  = "/home/tapas/ResultsOf"+TypeCal+"/"
```

2. In **support.py**

   - Update system dependent rotational B constant in **GetBconst()** function. It is needed only for linear rotor.

   - In **jobstring_sbatch()** function, adjust thread and walltime format.
   
```
thread         = Number of thread. 
walltime       = "40-00:00" 
```

   - In case of feynman, user must comment out the following #SBATCH command in the above mentioned function
   
```   
#SBATCH --account=rrg-pnroy
```

3. In **inputFile.py**

   - Make a list of beads in Getbeads() function. List of beads is defined by list_nb. Here basically same beades will be used for rotational and translational motions. If the user wish to use different set of beads, the user should consult with the developer.

   - Make three lists for step_trans, level, step in GetStepAndLevel() function. step_trans and step are the translational and rotational Monte Carlo step size. level is used in Monte Carlo bisection move for translational motion and it is integer in nature. Be careful, the function always needs the lists of step_trans, level, stepi, even if the user does not allow translation or rotational motions simultaneously. As for example, for the rotational motions only, the acceptance ration will be affected by the list of step (defined for rotational motion) only. Therefor, the user could fill up the step_trans, level lists by any real and integer numbers, respectively.



Now the script files are ready to submit your jobs. To know the command line arguments, just type the following command in terminal

python script_submission_analysis_MoRiBS.py -h

Examples of command line arguments to submit the jobs are given below:

        python script_submission_analysis_MoRiBS.py -d 1.0 -R 6.0 -N 2 -Block 100000 -Pass 100 --ROTMOVE tau submission PIMC H2O H2O 0.0333333333"

                -d       dipole moment value
                -R       Inter molecular distance
                -N       Number of rotors
                -Block   Number of Blocks
                -Pass    Number of Pass
                last value corresponds to the fixed beta value 1/T;

Read the outputs printed on the screen.

To analyze the output data -

        python script_submission_analysis_MoRiBS.py -d 1.0 -R 6.0 -N 2 -Block 100000 -Pass 100 --ROTMOVE --preskip 10000 tau analysis PIMC H2O H2O 0.0333333333"

Read the outputs printed on the screen. Final output files will be saved in directory final_results_path.


***N.B.: Don't hesitate to email to the developer if you face any problem to submit your jobs or analyze the output files by the scripts.***

Email: tapascuchem@gmail.com
