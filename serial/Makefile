release:
	gcc -o parsim parsim.c init_particles.c utils.c simulation.c cells.c -I.-O2 -fopenmp -lm

debug:
	gcc -o parsim parsim.c init_particles.c utils.c simulation.c cells.c -DDEBUG_MODE -I. -Wall -Wextra -g -fopenmp -lm

profile:
	gcc -o parsim parsim.c init_particles.c utils.c simulation.c cells.c -I. -g -O2 -fopenmp -lm

clean:
	rm -f parsim