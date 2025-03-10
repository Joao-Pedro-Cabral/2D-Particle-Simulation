release:
	gcc -o parsim-omp parsim-omp.c init_particles.c utils.c simulation.c cells.c -I.-O2 -fopenmp -lm -mavx2 -mfma

debug:
	gcc -o parsim-omp parsim-omp.c init_particles.c utils.c simulation.c cells.c -DDEBUG_MODE -I. -Wall -Wextra -g -fopenmp -lm -mavx2 -mfma

profile:
	gcc -o parsim-omp parsim-omp.c init_particles.c utils.c simulation.c cells.c -I. -g -O2 -fopenmp -lm -mavx2 -mfma

clean:
	rm -f parsim-omp