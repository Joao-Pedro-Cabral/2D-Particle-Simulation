
#include "particles_buffer.h"
#include <stdlib.h>

void particles_buffer_init(particles_buffer_t *buffer, long num_counters) {
    buffer->particles = NULL;
    buffer->counters = (long long *)malloc(num_counters * sizeof(long long));
    buffer->num_counters = num_counters;
    buffer->size = 0;
    buffer->total_num_particles = 0;
    buffer->capacity = 0;
    particles_buffer_clean_counters(buffer);
}

void particles_buffer_increment_counter(particles_buffer_t * buffer, long cell) {
    buffer->counters[cell] ++;
    buffer->size ++;
    if (buffer->size > buffer->capacity) {
        long long capacity = (buffer->capacity == 0) ? 1 : buffer->capacity * 2;
        particles_buffer_realloc(buffer, capacity);
    }
}

int particles_buffer_size(particles_buffer_t* buffer) {
    return buffer->size*sizeof(long long) + buffer->capacity*sizeof(particle_t) + sizeof(particles_buffer_t);
}

void particles_buffer_realloc(particles_buffer_t* buffer, long long capacity) {
    buffer->particles = (particle_t *)realloc(buffer->particles, capacity * sizeof(particle_t));
    buffer->capacity = capacity; 
}

void particles_buffer_clean_counters(particles_buffer_t *buffer) {
    for(long i = 0; i < buffer->num_counters; i ++) {
        buffer->counters[i] = 0;
    }
}

void particles_buffer_clean(particles_buffer_t *buffer) {
    free(buffer->particles);
    free(buffer->counters);
    buffer->particles = NULL;
    buffer->counters = NULL;
    buffer->num_counters = 0;
    buffer->size = 0;
    buffer->total_num_particles = 0;
    buffer->capacity = 0;
}
