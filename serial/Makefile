release:
	g++ -o parsim parsim.cpp init_particles.c utils.c simulation.cpp particles.cpp -I.-O2 -fopenmp

debug:
	g++ -o parsim parsim.cpp init_particles.c utils.c simulation.cpp particles.cpp -DDEBUG_MODE -I. -Wall -Wextra -g -fopenmp

clean:
	rm -f parsim