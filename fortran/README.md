# MC/MD Simulations

- To compile the code and generate MC/MD simulation data :

    ``` bash
        $ bash compile.sh
        # Generates common random initial coordinates and velocities
        $ bash gen_input.sh
        $ bash run_mc.sh &> mc_run.log
        $ bash run_md.sh &> md_run.log
    ```

- To calculate the radial distribution function from the `.xyz` file use `radialdist.x` as follows:

    ```sh
    $ radialdist.x -xyz md_lj.xyz -o md_lj.rdf -dr 0.02
    ```

    Summary of flags to use for command-line interface (cli):

    ```markdown
    - xyz    => .xyz positions trajectory filename
    - o      => .rdf output filename
    - dr     => bin-width
    - nstart => no. time-steps after which to start computation.
    ```

- To calculate time-correlation functions from the `.xyz` and `.vel` files use `diff.x` as follows:

    ```sh
     $ diff.x -xyz md_lj.xyz -vel md_lj.vel -o md_lj.diff -delta 0.1 -nt 50 -dt0 10
    ```

    Summary of flags to use for command-line interface (cli):

    ```markdown
    - xyz    => .xyz positions trajectory filename
    - vel    => .vel velocities trajectory filename
    - o      => .rdf output filename
    - delta  => time-difference between time-steps in traj. file.
    - nt     => Max correlation time (as a multiple of delta)
    - dt0    => time interval to store origins
    - nstart => no. time-steps after which to start computation.
    ```

## Unit Conventions

- Following are some useful quantities for conversion between real and reduced units, in case of Argon gas.

> - $\rho^{*} = 1.0 \leftrightarrow \rho = 1680 \text{ kg}/\text{m}^3$
> - $\text{T}^{*} = 1.0 \leftrightarrow \text{T} =  119.8 \text{ K} $
> - $\text{P}^{*} = 1.0 \leftrightarrow \text{P} =  41.9 \text{ MPa} $
> - $\Delta t^{*} = 1.0 \leftrightarrow \Delta t= 2.18 \text{ ps}$
> - $k_b = 1.38 \times 10^{-23} \text{ J}\text{K}^{-1}$
> - $\sigma = 3.405 \AA$
> - $\epsilon = 0.996 \text{ KJ }\text{mol}^{-1}$

- All coordinates in `.xyz` files are stored in angstroms ($\AA$) and velocities in `.vel` files are stored in angstroms per picoseconds ( $\AA/\text{ps}$ )

## File structure

- Root directory:   

```md
├── bin             => executable binaries
├── data            => MC/MD simulation data
├── images          => plot images from analysis
├── input_files     => sample input files
├── notebooks       => ipython notebooks with data analysis
├── src             => dir. with source code
└── report          => all report files
```

- Source code:  

```md
└── src
    ├── io_module.f90    => module to handle i/o
    ├── initialize.f90   => module for generating initial conditions
    ├── montecarlo.f90   => MC simulations module
    ├── moldyn.f90       => MD simulations module 
    ├── init.f90         => Generates positions and velocities
    ├── mc_nvt_lj.f90    => Runs NVT MC simulation
    ├── md_nve_lj.f90    => Runs NVE MD simulation
    ├── radialdist.f90   => Calculates radial distribution function
    └── diffusion.f90    => Calculates x,v correlation functions  
```

## Input Files

- Input parameters for generating initial co-ordinates and velocities can be passed to `init.x` through a standard namelist `&init_params` as follows:

```f90
&init_params

positions          = (logical) ! switch for generating initial positions

velocities         = (logical) ! switch for generating initial velocities

element            = (char=2)  ! type of atoms/particles.

n_part             = (integer) ! no. of particles

xyz_file           = (char=*)  ! file to store initial xyz configurations 

rho                = (float)   ! number density in the cell

sigma_angs         = (float)   ! min. distance b/w particles (in Angstroms)

fill_type          = (char=*)  ! type of initial config.: random, cubic, fcc, bcc

velocity_file      = (char=*)  ! file to store initial velocities

mass               = (float)   ! mass of the particles 

temperature        = (float)   ! temperature of initial particle

sample_type        = (char=*)  ! sampling v-distribution : uniform, boltzmann

/

```

- Input parameters for running Monte Carlo simulations can be passed to  `mc_nvt_lj.x` through a standard namelist `&mc_params` as follows:

```f90

&mc_params

initial_positions   = (char=*)  ! initial configuration file

outfile             = (char=*)  ! output file where energies are to be printed. 

xyz_file            = (char=*)  ! file to store xyz configurations 

temperature         = (float)   ! temperature of simulation 

delta_max           = (float)   ! maximum allowed displacement

ncycle              = (integer) ! total no. of Monte Carlo simulation steps

print_nstep         = (integer) ! prints xyz/observables every n-steps

sigma_angs          = (float)   ! min. distance b/w particles (in Angstroms) 

eps_kjpm            = (float)   ! epsilon energy parameter (in KJ/mol)

r_cut               = (float)   ! cut-off distance

/   

```

- Input parameters for running Molecular Dynamics simulations can be passed to  `md_nve_lj.x` through a standard namelist `&md_params` as follows:

```f90

&md_params

initial_positions   = (char=*)  ! initial configuration file

initial_velocities  = (char=*)  ! initial velocities file

outfile             = (char=*)  ! output file where energies are to be printed. 

xyz_file            = (char=*)  ! file to store xyz configurations 

velocity_file       = (char=*)  ! file to store velocities 

integration         = (char=*)  ! methods: verlet, velocity-verlet

ti                  = (float)   ! initial time 

tf                  = (float)   ! final time

del_t               = (float)   ! time step

print_nstep         = (integer) ! print observables for every n-steps

temperature         = (float)   ! temperature of simulation

sigma_angs          = (float)   ! min. distance b/w particles (in Angstroms) 

eps_kjpm            = (float)   ! epsilon energy parameter (in KJ/mol)

r_cut               = (float)   ! cut-off distance

/   

```


