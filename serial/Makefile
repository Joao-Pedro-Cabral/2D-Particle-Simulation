release:
	g++ -o parsim main.cpp init_particles.c utils.c parsim.cpp particles.cpp -I.-O2 -fopenmp

debug:
	g++ -o parsim main.cpp init_particles.c utils.c parsim.cpp particles.cpp -DDEBUG_H -I. -Wall -Wextra -g -fopenmp

clean:
	rm -f parsim