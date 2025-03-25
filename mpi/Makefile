release:
	mpicc -o parsim-mpi parsim-mpi.c init_particles.c utils.c simulation.c cells.c communication_buffers.c particles_buffer.c -I. -g -O2 -fopenmp -lm -mavx2 -mfma 

debug:
	mpicc -o parsim-mpi parsim-mpi.c init_particles.c utils.c simulation.c cells.c communication_buffers.c particles_buffer.c -DDEBUG_MODE -I. -Wall -Wextra -Wunused -g -fopenmp -lm -mavx2 -mfma 

profile:
	mpicc -o parsim-mpi parsim-mpi.c init_particles.c utils.c simulation.c cells.c communication_buffers.c particles_buffer.c -I. -g -O2 -fopenmp -lm -mavx2 -mfma 

clean:
	rm -f parsim-mpi