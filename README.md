Here are few instructions that will help a user to modify script files before submitting jobs in queue. The source codes can be downloaded directly from github by typing the following line on the terminal:

        git clone git@github.com:tapassahoo/MoRiBS-PIGS.git

First, the author must read README file in MoRiBS-PIGS/ and follow the instructions.

In the source directory, there are many Makefiles. Makefile-PIMC and Makefile-PIGS are the makefiles that a user needs to compile the source codes for finite temperature (PIMC) and ground state (PIGS) canculations, respectively. To compile the source codes, first copy Makefile-PIMC of Makefile-PIGS to Makefile and use the following command:
                               make clean
                               make

PotFunc() function in mc_estim.cc file includes analytic potential like dipole-dipole interaction.

But the user does not need to compile the source codes manually if the user like to submit jobs by script files. The script files are in dir: MoRiBS-PIGS/examples/scripts). There are three python scripts:

A. script_submission_analysis_MoRiBS.py
B. support.py
C. inputFile.py

The user are suggested to make the following modifications in the scripts before running MoRiBs successfully:

#------------------------------------------------------------------------#

A. In script_submission_analysis_MoRiBS.py

1. If user wish to run MoRiBs in graham.computecanada.ca, just replace "NameOfServer = "nlogn"" by "NameOfServer = "graham"". "NameOfServer = "nlogn"" when jobs will be submitted in feynman or nlogn server.

2. As computation of rotational density matrix for linear rotors is not time consuming, the user could use "status_rhomat = "Yes"" flag that generates rotational density matrix files instantly during submitting the jobs in queue. On the other hand, the user is advised to generate the density matrices for non linear molecules before he/she submit the jobs by the scripts. For this case, the user should use "status_rhomat = "No"".

3. If the user wish to include cage potential, he/she should use "status_cagepot = "Yes", otherwise, "status_cagepot = "No".

4. Keep the same directotory-tree as as the developer used - /home/user_name/source_dir/input_dir. The user may change the names of the directories. As for example, the developer used
        user_name           = "tapas"
        source_dir          = "Moribs-pigs/MoRiBS-PIMC/"                    #Path of the source directory#
        out_dir             = "nonlinear-molecule/"                         #This directory automatically generated in /work or in /scratch.
        input_dir           = "examples/nonlinear-molecule/"                #Where all the input and scripts are
        final_results_path  = "/home/"+user_name+"/ResultsOf"+TypeCal+"/"   #Where all the final results will be stored after analyzing the MoRiBs outputs.

but the user may change these as
        user_name           = "user_name" excluding "/"
        source_dir          = "MoRiBS-PIMC/"
        out_dir             = "PIMC-H2O/"
        input_dir           = "INPUT/"
