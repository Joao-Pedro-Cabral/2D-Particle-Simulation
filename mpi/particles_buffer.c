
#include "particles_buffer.h"
#include <stdlib.h>

void particles_buffer_init(particles_buffer_t *buffer, long num_counters, long long capacity) {
    buffer->particles = (particle_t *)malloc(capacity * sizeof(particle_t));
    buffer->counters = (long long *)malloc(num_counters * sizeof(long long));
    buffer->num_counters = num_counters;
    buffer->capacity = capacity;
    particles_buffer_clean_counters(buffer);
}

void particles_buffer_increment_counter(particles_buffer_t * buffer, long cell) {}

void particles_buffer_clean_counters(particles_buffer_t *buffer) {
    for(long i = 0; i < buffer->num_counters; i ++) {
        buffer->counters[i] = 0;
    }
}
void particles_buffer_add(particles_buffer_t* buffer, long cell) {}
void particles_buffer_clean(particles_buffer_t *buffer) {
    free(buffer->particles);
    free(buffer->counters);
    buffer->particles = NULL;
    buffer->counters = NULL;
    buffer->num_counters = 0;
    buffer->capacity = 0;
}
