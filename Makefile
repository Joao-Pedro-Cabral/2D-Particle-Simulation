release:
	gcc -o parsim main.c init_particles.c utils.c -I. -lm -O2

debug:
	gcc -o parsim main.c init_particles.c utils.c -DDEBUG_H -I. -lm -Wall -Wextra -g

clean:
	rm -f parsim