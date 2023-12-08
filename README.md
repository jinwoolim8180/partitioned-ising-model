# Fortran Implementation of Ising model in 2D

A Fortran implementation of Markov Chain Monte Carlo for 2D square-lattice Ising model.

It is possible to calculate mean energy, magnetization, specific heat, and susceptibility at various temperatures and save it to a csv file.

- [x] OpenMP supported. This repository is useful for simulations on a large square-lattice at few temperature points.

> [!Warning]
> For precise results, experiments on a large scale lattice Ising model need *a lot of* energy and time. We strongly recommend you to use a server with decent multi-core CPUs.

<!-- ## Result -->
<!-- 
![Result of 3D lattice](./result/plot_L30_D3_EQ16000_MC16000.png) -->

## How to run Monte Carlo simulation

### Clone the git repository

First, clone the git repository on your own computer directory (e.g. home directory).

After cloning the repository, enter the source directory.

```bash
git clone https://github.com/ising-model/ising-model.git
cd ising-model
```

### Compile the code

In order to run the simulation, you have to compile the code first.

To compile the code, simply run the command below:

```bash
bash compile.sh
```

### Run the simulation

To run experiments, run the command below:

```bash
./main --size 30 --init_temp 1.5 --final_temp 3.5 --temp_step 0.02 --eqstep 1000 --mcstep 1000
```

#### Options

To run simulation with your own custom options, run the program with the options below:

- -s, --block_size :  size of the partitioned block of a  lattice (default: 30)
- -i, --init_temp :   initial temperature of the output (default: 1.5)
- -f, --final_temp :  final temperature of the output   (default: 6.5)
- -t, --temp_step :   step size of the temperature      (default: 0.04)
- -m, --mcstep :      number of Monte Carlo steps       (default: 1000)
- -e, --eqstep :      number of steps for equilibration (default: 1000)
- -r, --dir :         directory to save the results     (default: ./results/)
- -p, --thread_per_row : number of threads per row of a lattice (default: 6)
- -h, --help :        print usage information and exit

## Future works to be done
We want to parallelize the sampling procedure using GPU.

We also want to speed up the process using various techniques (e.g. importance sampling). 

If you have abundant knowledge of those techniques, please contact us!
