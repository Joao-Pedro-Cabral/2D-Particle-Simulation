
# 2D Particle Simulation

A high-performance 2D particle simulation framework implemented in C. This project serves as a research and experimentation platform for parallel and distributed computing techniques applied to particle dynamics simulations (n-body problem), leveraging the Particle-in-Cell (PIC) approximation for computational efficiency.

## Overview

This simulation models the movement of particles under gravitational forces. Three distinct implementations allow users to explore various parallelization strategies:
- **Serial Implementation**: Provides a baseline, optimized for cache locality and SIMD acceleration (AVX2).
- **OpenMP Implementation**: Uses shared-memory parallelism with dynamic scheduling and per-cell locks to ensure thread-safe operations.
- **MPI Implementation**: Employs distributed-memory techniques following Foster’s Design Methodology (Partitioning, Communication, Agglomeration, Mapping) for scalable multi-node execution.

Detailed performance analyses and design decisions, including numerical precision challenges and load balancing strategies, are discussed in the reports on `doc/`.

## Features

- **Multiple Execution Modes**: Serial, OpenMP, and MPI showcase different implementations that explores parallelism at multiple levels.
- **Advanced Parallelization Strategies**:
  - **OpenMP**: Uses a single parallel region, loop parallelization with dynamic scheduling, and reduction clauses to manage unbalanced workload and avoid data races.
  - **MPI**: Implements communication buffers with multiple locks per neighbor, checkerboard-style domain agglomeration, and computation–communication overlap to reduce latency.
- **Precision and Performance Trade-offs**: The serial and OpenMP implementations use AVX2 SIMD intrinsics to boost performance despite introducing numerical precision due to the reordering of floating-point operations.
- **Testing and Benchmarking Suite**: Contains a dedicated `test-suite-mpi/` directory and sample Slurm job scripts for performance evaluation across various configurations.

## Implementation Details

### Serial Implementation

- **Algorithm Overview**: The simulation divides the 2D domain into a grid of cells. For each simulation step, the following are computed sequentially:
  - **Centers of Mass**: Each cell computes the weighted center from its particles.
  - **Kinetics**: The forces from particles in the same and neighboring cells are calculated to update positions and velocities.
  - **Particle Reassignment**: Particles moving out of their current cell are reallocated accordingly.
  - **Collision Handling**: Detects and processes particle collisions based on a predefined proximity threshold.
- **Optimization Techniques**: Utilizes a custom cell-based data structure that groups particle properties in vectors to improve cache usage and leverages AVX2 intrinsics for SIMD acceleration. Note that using SIMD reorders operations, which may lead to small discrepancies in floating-point arithmetic over many iterations.

### OpenMP Implementation

- **Parallel Design**:
  - **Single Parallel Region**: Minimizes overhead by spawning threads once outside of the timestep loop.
  - **Loop-Level Parallelism**: Uses `#pragma omp for` to distribute cell computations across threads.
  - **Dynamic Scheduling**: Balances workloads across cells with varying particle counts.
  - **Thread-Safety**: Implements per-cell locks to ensure safe concurrent updates while reducing contention.
  - **Reduction Clauses**: Utilizes OpenMP reductions for efficient collision count accumulation.
- **Performance Notes**: Benchmarks reveal significant speedups with 2, 4, and 8 threads.

### MPI Implementation

- **Design Strategy**: Follows Foster’s Design Methodology:
  - **Partitioning**: Divides the simulation grid into cells.
  - **Communication**: Implements dedicated buffers and locking mechanisms for safe particle exchange between neighboring grid partitions.
  - **Agglomeration**: Combines cells into larger blocks (via checkerboard or row-wise grouping) to reduce the communication overhead.
  - **Mapping**: Distributes cell blocks evenly across processors ensuring a balanced workload.
- **Optimizations**:
  - **Pre-Allocation of Buffers**: Reduces runtime overhead by avoiding repeated memory allocations and by initializing MPI requests.
  - **Communication–Computation Overlap**: Starts MPI receive requests early so computation can proceed concurrently.
- **Testing and Evaluation**: Test cases were conducted over varying numbers of MPI tasks and CPUs per task. Detailed results concerning execution time trends and speedup factors are presented in the project reports at `doc/`.

## Getting Started

### Prerequisites

- **C Compiler**: Ensure you have a C compiler installed (e.g., `gcc`).
- **AVX2**: Ensure your C compiler and your machine supports **AVX2**.
- **Make**: Required for building the project.
- **MPI Library**: For the MPI mode, install an MPI implementation like OpenMPI.
- **OpenMP Support**: Verify that your compiler supports OpenMP.
- **Slurm (Optional)**: For executing MPI jobs in a clustered environment using provided job scripts.

### Building the Project

Navigate to the desired execution mode directory and build the project:

#### Serial Version

```bash
cd serial
make
```

#### OpenMP Version

```bash
cd omp
make
```

#### MPI Version

```bash
cd mpi
make
```

### Running the Simulation

After building, run the simulation as follows:

#### Serial Version

```bash
./particle_simulation
```

#### OpenMP Version

```bash
./particle_simulation
```

#### MPI Version

```bash
mpirun -n <number_of_processes> ./particle_simulation
```

Replace `<number_of_processes>` with the desired number of MPI processes.

### Testing

A comprehensive test suite is provided in the `test-suite-mpi/` directory. For MPI testing on a cluster, refer to the included Slurm job submission scripts. An example script snippet:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=ParticleSim
#SBATCH --output=sim_output.txt
#SBATCH --error=sim_error.txt
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
srun ./particle_simulation [params]
```

## Project Structure

```
2D-Particle-Simulation/
├── serial/             # Serial implementation (AVX2 optimized)
├── omp/                # OpenMP parallel implementation
├── mpi/                # MPI distributed implementation
├── test-suite-mpi/     # Testing suite and job scripts for MPI configurations
├── .gitignore
└── README.md
```

## Performance and Evaluation

- **Execution Times & Scalability**: The report documents how different test cases scale with varying numbers of tasks and CPUs per task, noting near-linear scalability in the MPI version for up to 32 tasks.
- **Speedup Analyses**: Detailed measurements (see tables and figures in the reports) highlight the trade-offs between increased threading and overhead from synchronization, memory contention, and communication latency.
- **Known Issues**: Numerical instabilities from AVX2-based SIMD optimizations are acknowledged—minor floating-point discrepancies may result in slight variations in collision counts between runs.

For further information on benchmarks, detailed methodology, and performance graphs, please refer to the attached reports on `doc/`.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with any enhancements or bug fixes. When contributing, include updates to the documentation or tests as needed.

## Future Enhancements

- **Dynamic Load Balancing**: Implement strategies to further balance workloads between nodes dynamically.
- **Enhanced Precision Handling**: Explore alternative algorithms or compensations to mitigate the numerical precision issues inherent in SIMD operations.
- **Extended Benchmarking**: Automate and expand the performance testing across diverse hardware configurations.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

Developed by João Pedro Cabral, Lisa Santarossa and Matilde de Almeida e Paiva Vital Ferreira. For additional projects and updates, visit [João Pedro Cabral's GitHub](https://github.com/Joao-Pedro-Cabral).
