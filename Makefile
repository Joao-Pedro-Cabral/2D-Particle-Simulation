release:
	g++ -o parsim main.cpp init_particles.c utils.c parsim.cpp particles.cpp -I. -lm -O2 -fopenmp

debug:
	g++ -o parsim main.cpp init_particles.c utils.c parsim.cpp particles.cpp -DDEBUG_H -I. -lm -Wall -Wextra -g -fopenmp

clean:
	rm -f parsim