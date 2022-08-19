## User's Guide

**Note:** It is derived from the MoRiBS (see `Computer Physics Communications vol. 204, pp. 170-188, 2016`) to estimate ground state properties of many-body interacting quantum systems based on the Path Integral Ground State (PIGS) approach.

**Disclaimer**: We disclaim all warranties concerning the enclosed program.

The source codes of `MoRiBS-PIGS` are mainly written by Tapas Sahoo and Pierre-Nicholas Roy, University of Waterloo, Canada. If you encounter any problem with compiling or running the program, don't hesitate to contact Tapas Sahoo by email `tapascuchem@gmail.com`.

### **Setup**

To download the package, a user should install `git` on his/her computer, followed by running the below command on his/her terminal:
```
git clone https://github.com/tapassahoo/MoRiBS-PIGS.git
```
***Users should not change the directory structure after they unzip the distributed file.*** We label the main directory, where the source codes (`*.cc`, `*.h`, and `*.f`) are, as `$MAIN`. Please note that we have provided example configuration files like `Makefile` with the distribution and one may just adjust the files according to the architecture of his/her computer. There is no need to create new configuration files. 

### **Compilation** 

The compilation procedure of `MoRiBS-PIGS` contains the following steps:

1. `cd $MAIN/spring` and open `make.CHOICES`. One should specify the platform of his/her computer by uncommenting the correct line, i.e., `PLAT = LINUX`;
2. `cd $MAIN/spring/SRC` and open `make.${PLAT}`. With the above choice, it should be `make.LINUX`. In `make.${PLAT}`, one should specify the `fortran` and `C/C++` compilers that are installed in his/her computer. Note that only the non-MPI compilation of SPRNG has been tested with MoRiBS-PIGS;
3. at `$MAIN/spring/SRC`, `make clean` and `make`. The generated libraries are in `$MAIN/spring/lib`. Make sure `$MAIN/spring/lib` is empty before `make`. Sometimes, `make clean` may not clean up `$MAIN/spring/lib` completely. Steps 1-3 are to compile the sprng libraries that are needed by the main code;
4. at `$MAIN`, open `Makefile`, and specify options, `CFLAGS`, `LDFLAGS`, `CC` and `FC`. Those are the optimization options, `C/C++` compilation flags, link flags, `C/C++` compiler, and `Fortran` compiler respectively;
5. at `$MAIN`, `make -f Makefile-GNU clean` and `make -f Makefile-GNU`. If there is no error message, the generated executable is called `pimc`.


					************************
					**   BEFORE RUNNING   **
					************************

Please note that one should delete the following files before each simulation run:
	1: yw001.*;
	2: permutation,
Otherwise, the simulation will bomb out.

Remember to use asymrho.x (in nmv_prop/), symrho.x (in symtop_prop/) or linden.x (in linear_prop/) to generate the needed files for PIMC sampling of asymmetric top, symmetric top, or linear rotors. Users should read the respective README files in those directories carefully before compiling and running those programs. Users should carefully name the resultant files in accordance to the rules explained in Computer Physics Communications vol. 204, pp. 170–188, yr 2016.

N.B.: For the MoRiBS-PIGS simulations of the nonlinear rotors, the user should always use an odd number of beads, and the name of all *.rho, *.eng, *.esq files must be of the format - {rotor name}_T{temperature}t{# of rotational beads}.rho/eng/esq. Remember, the {temperature} string's total width must be of six digits after the decimal point.

As for Example:

	H2O_T10.000000t101.rho	
	H2O_T10.000000t101.eng
	H2O_T10.000000t101.esq

To generate the above files by asymrho.f, P-1 beads are used in MoRiBS-PIGS simulations instead of P beads used in the MoRiBS-PIMC simulations. Imaginary time (tau) is defined as (beta/P-1), where P is the number of beads.

The name of the working directory to submit jobs is ~/MoRiBS-PIGS/examples/pigs-simulations-for-para-Water. To execute the simulations, the user must have the following files -
	1: qmc.input
	2: *.rho, *.eng, *.esq files

The user is advised to read Computer Physics Communications vol. 204, pp. 170–188, yr 2016, to set up the input (qmc.input).

One may compile the code, move the resultant executable "pimc" to examples/MF_8He_0.37K_512_128 or examples/N2O_5pH2_512_128, and try running the respective simulations. One may also look at the files in those two example directories to have a sense of how many files are needed and their formats. Some of these examples are discussed in our associated paper.

					************************
					**   OUTPUT SUMMARY   **
					************************

There are two relevant directories for a MoRiBS-PIGS run, the working directory and the output directory. The former is the directory where the program executable is and the latter is specified in qmc.input. In below we simply call the two directories $WORK and $OUTPUT. And we also use $FILNAM to specify the FILENAMEPREFIX that is specified in qmc.input.

************
   energy
************
All energy outputs are stored in the following two files:

	$OUTPUT/$FILNAM_sum.eng and $OUTPUT/$FILNAM.eng.

The first contains the energy components averaged on the fly while the latter contains the block-averaged quantities for every block. In $OUTPUT/$FILNAM_sum.eng, column
	1: counting indices of punching out the energy components;
	2: potential energy;
	3: total energy (rotational + potential);

In $OUTPUT/$FILNAM.eng, column
	1: block index;
	2: block averaged potential energy;
	3: block averaged total energy (rotational + potential);

All energy quantities have the unit of K.

*******************
   distribution
*******************

Another important output is output.xyz where all the instantaneous values of rotational degrees of freedom in the space-fixed frame are stored. For better understanding, the user can see the SaveInstantConfig() function in mc_main.cc.

	loop over rotors' index
		loop over rotational beads
			loop over (cos(theta),phi,chi)

The output.xyz file is required to plot order parameters, angular distribution functions and <O-H..O distributions.
 

