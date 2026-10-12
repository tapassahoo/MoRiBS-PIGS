# User's Guide to MoRiBS-PIGS

## 1. Overview

MoRiBS-PIGS is a computational package derived from MoRiBS (Computer
Physics Communications, **204**, 170--188, 2016) for estimating
ground-state properties of many-body interacting quantum systems using
the Path Integral Ground State (PIGS) method.

The source code of MoRiBS-PIGS was developed primarily by Tapas Sahoo
and Pierre-Nicholas Roy at the University of Waterloo, Canada.

**Disclaimer:** The authors disclaim all warranties concerning the
enclosed software.

For assistance with compilation or execution, please contact Tapas Sahoo
at **tapascuchem@gmail.com**.

## 2. Installation and Setup

### 2.1. Downloading the Source Code

Before downloading the package, ensure that Git is installed on your
computer. Open a terminal and execute:

``` bash
git clone https://github.com/tapassahoo/MoRiBS-PIGS.git
```

**Important:** Preserve the original directory structure of the
downloaded package. Do not move or rename directories unless explicitly
required by the configuration files.

Throughout this guide, `$MAIN` denotes the main directory of the
MoRiBS-PIGS package, which contains the principal source files (`*.cc`,
`*.h`, and `*.f`).

The distribution includes example configuration files, such as
`Makefile`, which can be modified to suit the user's computer
architecture and installed compilers. Creating new configuration files
is generally unnecessary.

## 3. Compilation

The compilation procedure consists of the following steps.

### Step 1: Select the Platform

Navigate to the `spring` directory:

``` bash
cd $MAIN/spring
```

Open `make.CHOICES` and uncomment the line corresponding to your
operating system and platform. For example, for Linux, select:

``` makefile
PLAT = LINUX
```

Ensure that the selected platform is compatible with your system.

### Step 2: Configure the SPRNG Compilation

Navigate to the SPRNG source directory:

``` bash
cd $MAIN/spring/SRC
```

Open the platform-specific configuration file. For Linux, this file is
`make.LINUX`.

Specify the Fortran and C/C++ compilers available on your system,
together with any necessary compiler options.

**Note:** Only the non-MPI version of SPRNG has been tested with
MoRiBS-PIGS.

### Step 3: Compile the SPRNG Libraries

Before compilation, ensure that `$MAIN/spring/lib` is empty. Then
execute:

``` bash
make clean
make
```

The compiled SPRNG libraries should be generated in:

``` text
$MAIN/spring/lib
```

Check that the directory contains the expected libraries before
proceeding. In some cases, `make clean` may not remove all previously
generated files from this directory; therefore, verify its contents
manually.

Steps 1--3 prepare the SPRNG libraries required by the main program.

### Step 4: Configure the Main Program

Navigate to the main directory:

``` bash
cd $MAIN
```

Open `Makefile` and configure the following variables according to your
system:

-   `CFLAGS`: Compiler options for C/C++ source files, including
    optimization and other compilation flags.
-   `LDFLAGS`: Linker options and required library paths.
-   `CC`: C/C++ compiler.
-   `FC`: Fortran compiler.

Ensure that the selected compilers and flags are compatible with the
source code and the compiled SPRNG libraries.

### Step 5: Compile MoRiBS-PIGS

From the main directory, execute:

``` bash
make -f Makefile-GNU clean
make -f Makefile-GNU
```

If compilation completes successfully, the executable named `pimc`
should be generated in the main directory.

## 4. Preparing and Running a Simulation

### 4.1. Cleaning Up Previous Simulation Files

Before starting each new simulation, remove the following files from the
working directory, as applicable:

``` bash
rm -f yw001.*
rm -f permutation
```

These files may contain data from a previous run. Leaving them in place
can interfere with a subsequent simulation.

**Caution:** Execute these commands only in the intended simulation
working directory. Verify the directory and the files before deleting
them.

### 4.2. Preparing Rotational Density-Matrix Files

Use the appropriate utility to generate the files required for PIGS
sampling, depending on the molecular rotor type:

  Rotor type       Utility       Directory
  ---------------- ------------- ----------------
  Asymmetric top   `asymrho.x`   `nmv_prop/`
  Symmetric top    `symrho.x`    `symtop_prop/`
  Linear rotor     `linden.x`    `linear_prop/`

Read the corresponding `README` file in each directory carefully before
compiling or executing these utilities.

The generated files must follow the naming conventions and formatting
requirements described in the associated MoRiBS publication (Computer
Physics Communications, **204**, 170--188, 2016).

### 4.3. Naming Conventions for Nonlinear Rotors

For MoRiBS-PIGS simulations involving nonlinear rotors, the number of
imaginary-time beads must be **odd**.

The rotational density-matrix and energy files must follow this naming
convention:

``` text
{rotor_name}_T{temperature}t{number_of_beads}.{extension}
```

Here, `{extension}` is `rho`, `eng`, or `esq`, and the temperature must
be represented using exactly six digits after the decimal point.

For example, for an H2O rotor at 10 K with 101 beads, the files should
be named:

``` text
H2O_T10.000000t101.rho
H2O_T10.000000t101.eng
H2O_T10.000000t101.esq
```

### 4.4. Imaginary-Time Discretization

For generating the files used in MoRiBS-PIGS simulations with
`asymrho.f`, use (P-1) beads in the density-matrix preparation, rather
than the (P) beads used in the corresponding MoRiBS-PIMC procedure.

The imaginary-time step is defined as

\[ `\tau`{=tex}=`\frac{\beta}{P-1}`{=tex}, \]

where (`\beta`{=tex}) is the inverse temperature in the convention used
by the program and (P) is the number of beads in the PIGS simulation.

**Note:** The bead-count convention and the corresponding imaginary-time
discretization must be consistent with the implementation described in
the associated publication.

### 4.5. Required Input Files

For the example simulation in

``` text
~/MoRiBS-PIGS/examples/pigs-simulations-for-para-Water
```

the working directory must contain, at a minimum:

1.  `qmc.input`
2.  The required rotational density-matrix and energy files (`*.rho`,
    `*.eng`, and `*.esq`).

The input file `qmc.input` specifies the simulation parameters. Consult
the associated publication for guidance on selecting these parameters
and preparing a valid input file.

For additional guidance, examine the example directories:

``` text
examples/MF_8He_0.37K_512_128
examples/N2O_5pH2_512_128
```

After successful compilation, the `pimc` executable can be copied into
the appropriate example directory and executed there. These examples
illustrate the required input files, directory organization, and file
formats; some are also discussed in the associated publication.

## 5. Output Files and Their Interpretation

A MoRiBS-PIGS simulation involves two principal directories:

-   **Working directory (`$WORK`):** The directory from which the
    executable is run and in which the simulation is configured.
-   **Output directory (`$OUTPUT`):** The directory specified in
    `qmc.input` for storing simulation results.

The variable `$FILNAM` denotes the filename prefix specified by
`FILENAMEPREFIX` in `qmc.input`.

### 5.1. Energy Outputs

The simulation produces two principal energy-output files:

``` text
$OUTPUT/$FILNAM_sum.eng
$OUTPUT/$FILNAM.eng
```

**File 1: `$FILNAM_sum.eng`**

This file contains energy components averaged during the simulation. Its
columns are interpreted as follows:

  Column   Quantity
  -------- ----------------------------------------------------------------
  1        Index associated with the recorded energy output
  2        Potential energy
  3        Total energy, including rotational and potential contributions

**File 2: `$FILNAM.eng`**

This file contains block-averaged energy quantities:

  -----------------------------------------------------------------------
  Column                              Quantity
  ----------------------------------- -----------------------------------
  1                                   Block index

  2                                   Block-averaged potential energy

  3                                   Block-averaged total energy,
                                      including rotational and potential
                                      contributions
  -----------------------------------------------------------------------

All energy quantities are expressed in **kelvin (K)**.

### 5.2. Rotational Configurations and Angular Distributions

Another important output file is:

``` text
output.xyz
```

This file stores instantaneous configurations of the rotational degrees
of freedom in the space-fixed frame.

Conceptually, the data are organized by looping over the rotor index,
then over the rotational beads, and finally over the angular
coordinates:

\[ (`\cos`{=tex}`\theta`{=tex},`\phi`{=tex},`\chi`{=tex}). \]

The precise output layout should be checked against the implementation
of the `SaveInstantConfig()` function in `mc_main.cc`.

The `output.xyz` file can be used to analyse rotational configurations
and calculate quantities such as:

-   Molecular orientational order parameters.
-   Angular distribution functions.
-   O--H···O distance distributions and related structural correlations,
    where applicable.

Appropriate post-processing is required to obtain these quantities from
the recorded configurations.

## 6. An Alternative Simulation Workflow

MoRiBS-PIGS can also be used through a Python-based execution workflow,
allowing users to run simulations without manually performing every
preparation and execution step.

The script

``` text
execution_driver_linear_molecule.py
```

provides an alternative entry point for simulations involving linear
molecules.

Before using the script, examine its command-line arguments,
configuration requirements, input-file conventions, and execution
instructions. Ensure that the required executable and input data are
available in the expected directories.

The precise workflow and supported options should be determined from the
script and its accompanying documentation.

## 7. Further Information

Users are strongly encouraged to consult the following publication for
the theoretical background, computational methodology, and detailed
simulation setup:

*Computer Physics Communications*, **204**, 170--188 (2016).

The example directories and the `README` files supplied with the package
provide additional practical guidance for preparing input files,
compiling auxiliary utilities, and executing simulations.

For questions concerning compilation, execution, or the use of
MoRiBS-PIGS, contact Tapas Sahoo at **tapascuchem@gmail.com**.
