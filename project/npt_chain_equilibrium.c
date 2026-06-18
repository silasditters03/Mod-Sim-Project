#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"
#include <stdlib.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 500

/* Initialization variables */
/* The simulation no longer targets a fixed number of MC steps.
   It runs until the volume has equilibrated.
   max_safety_steps is only a protection against an infinite run. */
const int max_safety_steps = 5000000;
const int output_steps = 100;
const int progress_print_steps = 5000;
const double packing_fraction = 0.1; // Starting density
const double diameter = 1.0;
const double delta  = 0.1;
/* Volume change -deltaV, delta V */
const double deltaV = 2.0;

/* REMOVED const so we can loop over pressure */
double betaP;
const char* init_filename = "xyz_500.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];
double volume_list[900] = {}; 
double chain_len = 5.0; //<--------------------------- nieuwe var

/* Equilibration detection settings.
   We measure the volume every output_steps and compare two adjacent
   volume blocks. 50 samples = 50 * output_steps MC steps. */
const int eq_min_steps = 10000;          // do not even test before this many MC steps
const int eq_block_samples = 50;         // volume samples per block
const int eq_required_stable = 3;        // require this many stable checks in a row
const double eq_rel_tol = 0.001;         // 0.5 percent relative block-mean tolerance


//Self made variables for easier computing
double diameter_squared  = diameter*diameter;

/* Cell List Arrays */
int *head = NULL; // Size: total_cells
int *lscl = NULL; // Size: N (n_particles)
int n_cells[NDIM];
double cell_size[NDIM];
int total_cells;

// Call this in main() after read_data()
void allocate_cell_list(void) {
    lscl = (int*)malloc(N * sizeof(int));
}

//Generates random double between any max and min value
double get_random_double(double min, double max) {
    if (min > max) {
        double temp = min;
        min = max;
        max = temp;
    }
    return min + (dsfmt_genrand() * (max - min));
}

//Applies periodic boundary conditions on a coordinate
double apply_pbc(double position, double box_length) {
    if (position >= box_length) position -= box_length;
    if (position < 0.0) position += box_length;
    return position;
}

int get_cell_index(double rx, double ry, double rz) {
    int cx = (int)(rx / cell_size[0]);
    int cy = (int)(ry / cell_size[1]);
    int cz = (NDIM == 3) ? (int)(rz / cell_size[2]) : 0;

    // Safety checks for boundary float precision
    if(cx >= n_cells[0]) cx = n_cells[0] - 1;
    if(cy >= n_cells[1]) cy = n_cells[1] - 1;
    if(NDIM == 3 && cz >= n_cells[2]) cz = n_cells[2] - 1;

    int cell_idx = cx + cy * n_cells[0];
    if (NDIM == 3) cell_idx += cz * n_cells[0] * n_cells[1];
    
    return cell_idx;
}

void build_cell_list(void) {
    total_cells = 1;
    for (int d = 0; d < NDIM; d++) {
        n_cells[d] = (int)(box[d] / diameter);
        if (n_cells[d] < 3) n_cells[d] = 3; 
        cell_size[d] = box[d] / n_cells[d];
        total_cells *= n_cells[d];
    }

    // Dynamic allocation/reallocation of head array
    head = (int*)realloc(head, total_cells * sizeof(int));
    for (int c = 0; c < total_cells; c++) head[c] = -1;

    // Construct linked list using the helper function
    for (int n = 0; n < n_particles; n++) {
        int cell_idx = get_cell_index(r[n][0], r[n][1], (NDIM == 3) ? r[n][2] : 0.0);
        lscl[n] = head[cell_idx];
        head[cell_idx] = n;
    }
}

void update_cell_list_local(int target_particle, int old_cell, int new_cell) {
    if (old_cell == new_cell) return; 

    // --- 1. REMOVE FROM OLD CELL ---
    if (head[old_cell] == target_particle) {
        head[old_cell] = lscl[target_particle];
    } else {
        int current = head[old_cell];
        while (current != -1 && lscl[current] != target_particle) {
            current = lscl[current];
        }
        if (current != -1) {
            lscl[current] = lscl[target_particle];
        }
    }

    // --- 2. ADD TO NEW CELL ---
    lscl[target_particle] = head[new_cell]; 
    head[new_cell] = target_particle;       
}

/* Functions */
int change_volume(void){
    double dV = get_random_double(-1*deltaV, deltaV);
    
    double volume_old = 1.0;
    double box_old[NDIM];
    for (int i = 0; i < NDIM; i++){
        volume_old *= box[i];
        box_old[i] = box[i];
    }
    
    double volume_new = volume_old + dV;
    if (volume_new <= 0) return 0;

    double half_length[NDIM];
    for (int i = 0; i < NDIM; i++){
        box[i] = cbrt(volume_new);
        half_length[i] = 0.5 * box[i];
    }
    
    double scale_x = box[0] / box_old[0];
    double scale_y = box[1] / box_old[1];
    double scale_z = box[2] / box_old[2];
    
    for (int n = 0; n < n_particles; n++){
        r[n][0] *= scale_x;
        r[n][1] *= scale_y;
        r[n][2] *= scale_z;
    }
    
    if (volume_new < volume_old){ 
        for (int n = 0; n < n_particles; n++){
            double rx = r[n][0];
            double ry = r[n][1];
            double rz = r[n][2];
            for (int m = n + 1; m < n_particles; m++){
                double dx = rx - r[m][0];
                double dy = ry - r[m][1];
                double dz = rz - r[m][2];
                
                if (dx > half_length[0]) dx -= box[0];
                if (dx <= -half_length[0]) dx += box[0];
                
                double dx2 = dx * dx; 
                if (dx2 >= diameter_squared) continue;

                if (dy > half_length[1]) dy -= box[1];
                if (dy <= -half_length[1]) dy += box[1];
                
                double dy2 = dy * dy; 
                if (dx2 + dy2 >= diameter_squared) continue;

                if (dz > half_length[2]) dz -= box[2];
                if (dz <= -half_length[2]) dz += box[2];
                
                double dz2 = dz * dz;
                
                if (dx2 + dy2 + dz2 < diameter_squared) {
                    for (int p = 0; p < n_particles; p++){
                        for (int i = 0; i < NDIM; i++){
                            r[p][i] = r[p][i] * (box_old[i] / box[i]);
                        }
                    }
                    for (int i = 0; i < NDIM; i++) box[i] = box_old[i];
                    return 0;
                }
            }
        }
    }
    
    double acc_arg = -(betaP / (diameter_squared * diameter)) * (volume_new - volume_old) + n_particles * log(volume_new / volume_old);
    double acc_exp = exp(acc_arg);
    double acceptance = fmin(1.0, acc_exp);
    
    if(dsfmt_genrand() < acceptance){
        build_cell_list();
        return 1; 
    }
    else {
        for (int n = 0; n < n_particles; n++){
            r[n][0] *= 1/scale_x;
            r[n][1] *= 1/scale_y;
            r[n][2] *= 1/scale_z;
        }
        for (int i = 0; i < NDIM; i++) box[i] = box_old[i];
        return 0;
    }
}

void read_data(void){
    FILE *file = fopen(init_filename, "r");
    if (file == NULL){
        printf("File not found\n");
        exit(1);
    }

    if(fscanf(file, "%d", &n_particles) != 1){
        printf("%d\n",n_particles);
        fclose(file);
        exit(1);
    }

    for(int d = 0; d < NDIM; d++){
        double min_val, max_val;
        if (fscanf(file, "%lf %lf", &min_val, &max_val) == 2) {
            box[d] = max_val - min_val;
        }
        else {
            printf("Error: Failed to read box dimensions on line %d.\n", d + 2);
            fclose(file);
            exit(1);
        }
    }

    double temp_diameter;
    if (NDIM == 2){
        for(int n = 0; n < N; n++){
            if (fscanf(file, "%lf %lf %*lf %lf", &r[n][0], &r[n][1], &temp_diameter) != 3) {
                printf("Error: Failed to read coordinates for particle %d.\n", n);
                fclose(file);
                exit(1);
            }
        }
    }
    else if (NDIM == 3){
        for(int n = 0; n < N; n++){
            if (fscanf(file, "%lf %lf %lf %lf", &r[n][0], &r[n][1], &r[n][2], &temp_diameter) != 4) {
                printf("Error: Failed to read coordinates for particle %d.\n", n);
                fclose(file);
                exit(1);
            }
        }
    }
    else {
        printf("Error: Unsupported NDIM value.\n");
        fclose(file);
        exit(1);
    }

    fclose(file);
    printf("Successfully read %d particles from %s.\n", n_particles, init_filename);
}

int check_cell_overlaps(int target_particle, double rx, double ry, double rz, double half_length[NDIM]) {
    int cx = (int)(rx / cell_size[0]);
    int cy = (int)(ry / cell_size[1]);
    int cz = (NDIM == 3) ? (int)(rz / cell_size[2]) : 0;

    for (int ox = -1; ox <= 1; ox++) {
        int n_cx = (cx + ox + n_cells[0]) % n_cells[0];
        
        for (int oy = -1; oy <= 1; oy++) {
            int n_cy = (cy + oy + n_cells[1]) % n_cells[1];
            
            int z_loops = (NDIM == 3) ? 1 : 0;
            for (int oz = -z_loops; oz <= z_loops; oz++) {
                int n_cz = (NDIM == 3) ? (cz + oz + n_cells[2]) % n_cells[2] : 0;

                int neighbor_cell = n_cx + n_cy * n_cells[0];
                if (NDIM == 3) neighbor_cell += n_cz * n_cells[0] * n_cells[1];

                int m = head[neighbor_cell]; 
                
                while (m != -1) {
                    if (m != target_particle) {
                        double dx = rx - r[m][0];
                        if (dx > half_length[0]) dx -= box[0];
                        if (dx <= -half_length[0]) dx += box[0];
                        if (dx * dx >= diameter_squared) { m = lscl[m]; continue; }

                        double dy = ry - r[m][1];
                        if (dy > half_length[1]) dy -= box[1];
                        if (dy <= -half_length[1]) dy += box[1];
                        if (dx * dx + dy * dy >= diameter_squared) { m = lscl[m]; continue; }

                        double dz = 0.0;
                        if (NDIM == 3) {
                            dz = rz - r[m][2];
                            if (dz > half_length[2]) dz -= box[2];
                            if (dz <= -half_length[2]) dz += box[2];
                        }

                        if (dx * dx + dy * dy + dz * dz < diameter_squared) return 1; 
                    }
                    m = lscl[m]; 
                }
            }
        }
    }
    return 0; 
}


int move_particle_classic(void){
    int particle_index = (int) get_random_double(0, n_particles);
    
    double old_r[NDIM];
    for (int i = 0; i < NDIM; i++) old_r[i] = r[particle_index][i];
    int old_cell = get_cell_index(old_r[0], old_r[1], (NDIM == 3) ? old_r[2] : 0.0);

    for (int i = 0; i < NDIM; i++){
        r[particle_index][i] = apply_pbc(r[particle_index][i] + get_random_double(-delta, delta), box[i]); 
    }
    
    double rx = r[particle_index][0];
    double ry = r[particle_index][1];
    double rz = (NDIM == 3) ? r[particle_index][2] : 0.0;
    
    double half_length[NDIM];
    for (int i = 0; i < NDIM; i++) half_length[i] = 0.5 * box[i];

    if (check_cell_overlaps(particle_index, rx, ry, rz, half_length)) {
        for(int i = 0; i < NDIM; i++) r[particle_index][i] = old_r[i];
        return 0; 
    }

    int new_cell = get_cell_index(rx, ry, rz);
    update_cell_list_local(particle_index, old_cell, new_cell);
    
    return 1;
}


int move_particle(int particle, double x, double y, double z, //coords
    double deltax, double deltay, double deltaz){

    int new_particle;
    bool overlap(){
        for(int i = 0; i < n_particles; i++){
            if(i == particle){continue;}
            double dx = x - r[i][0];
            double dy = y - r[i][1];
            double dz = z - r[i][2];

            if (dx >  0.5*box[0]) dx -= box[0];
            if (dx < -0.5*box[0]) dx += box[0];
            if (dy >  0.5*box[1]) dy -= box[1];
            if (dy < -0.5*box[1]) dy += box[1];
            if (dz >  0.5*box[2]) dz -= box[2];
            if (dz < -0.5*box[2]) dz += box[2];

            double dist = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
            
            if(dist < diameter){
                new_particle = i;
                return true;
            }
        }
        return false;
    }

    x += deltax;
    y += deltay;
    z += deltaz;

    if(x > box[0]){x -= box[0];}
    if(x < 0.){x += box[0];}
    if(y > box[1]){y -= box[1];}
    if(y < 0.){y += box[1];}
    if(z > box[2]){z -= box[2];}
    if(z < 0.){z += box[2];}


    if(!overlap()){
        r[particle][0] = x;
        r[particle][1] = y;
        r[particle][2] = z;
        return particle;
    }

    return new_particle;
}

static inline double min_image(double dx, double L){
    if(dx >  0.5*L) dx -= L;
    if(dx < -0.5*L) dx += L;
    return dx;
}

static inline void wrap_coordinate(double *x, double L){
    if(*x >= L || *x < 0.0){
        *x -= L * floor(*x / L);
        if(*x >= L) *x -= L;
        if(*x < 0.0) *x += L;
    }
}

int move_particle_chain(void){

    const double sigma = diameter;
    const double sigma2 = sigma * sigma;
    const double eps = 1e-12;

    double remaining = chain_len;

    /* Pick random starting particle */
    int current_particle = (int)(n_particles * dsfmt_genrand());
    if(current_particle >= n_particles) current_particle = n_particles - 1;

    /* Pick random direction and normalize it */
    double ex, ey, ez, norm2, invnorm;

    do{
        ex = 2.0 * dsfmt_genrand() - 1.0;
        ey = 2.0 * dsfmt_genrand() - 1.0;
        ez = 2.0 * dsfmt_genrand() - 1.0;

        norm2 = ex*ex + ey*ey + ez*ez;
    } while(norm2 < eps);

    invnorm = 1.0 / sqrt(norm2);
    ex *= invnorm;
    ey *= invnorm;
    ez *= invnorm;

    while(remaining > eps){

        double best_s = remaining;
        int hit_particle = -1;

        double x = r[current_particle][0];
        double y = r[current_particle][1];
        double z = r[current_particle][2];

        int old_cell = get_cell_index(x, y, z);

        /*
           Use the cell list to search only particles that can possibly collide.
           A particle can only be hit during the remaining chain length if its
           current center is within sigma + remaining of the moving particle.
           Because cell_size is about sigma, this usually means only a few cells.
        */
        double search_radius = sigma + remaining;
        double search_radius2 = search_radius * search_radius;

        int cx = (int)(x / cell_size[0]);
        int cy = (int)(y / cell_size[1]);
        int cz = (NDIM == 3) ? (int)(z / cell_size[2]) : 0;

        if(cx >= n_cells[0]) cx = n_cells[0] - 1;
        if(cy >= n_cells[1]) cy = n_cells[1] - 1;
        if(NDIM == 3 && cz >= n_cells[2]) cz = n_cells[2] - 1;

        int ox_max = (int)ceil(search_radius / cell_size[0]);
        int oy_max = (int)ceil(search_radius / cell_size[1]);
        int oz_max = (NDIM == 3) ? (int)ceil(search_radius / cell_size[2]) : 0;

        /* Find nearest collision in the chain direction using cell-list candidates */
        for(int ox = -ox_max; ox <= ox_max; ox++){
            int ncx = (cx + ox + n_cells[0]) % n_cells[0];

            for(int oy = -oy_max; oy <= oy_max; oy++){
                int ncy = (cy + oy + n_cells[1]) % n_cells[1];

                for(int oz = -oz_max; oz <= oz_max; oz++){
                    int ncz = (NDIM == 3) ? (cz + oz + n_cells[2]) % n_cells[2] : 0;

                    int cell_idx = ncx + ncy * n_cells[0];
                    if(NDIM == 3) cell_idx += ncz * n_cells[0] * n_cells[1];

                    int j = head[cell_idx];
                    while(j != -1){

                        if(j != current_particle){
                            double dx = x - r[j][0];
                            double dy = y - r[j][1];
                            double dz = z - r[j][2];

                            dx = min_image(dx, box[0]);
                            dy = min_image(dy, box[1]);
                            dz = min_image(dz, box[2]);

                            double r2 = dx*dx + dy*dy + dz*dz;

                            /* Too far away to be reached during the remaining chain */
                            if(r2 <= search_radius2){

                                /*
                                   Collision condition:
                                   |rij + s e|^2 = sigma^2

                                   b = rij . e
                                   c = rij^2 - sigma^2
                                   s = -b - sqrt(b^2 - c)
                                */
                                double b = dx*ex + dy*ey + dz*ez;

                                /* If b >= 0, particle j is not in front of the moving particle */
                                if(b < 0.0){
                                    double c = r2 - sigma2;
                                    double discriminant = b*b - c;

                                    if(discriminant > 0.0){
                                        double s = -b - sqrt(discriminant);

                                        if(s > eps && s < best_s){
                                            best_s = s;
                                            hit_particle = j;
                                        }
                                    }
                                }
                            }
                        }

                        j = lscl[j];
                    }
                }
            }
        }

        /* Move current particle either to next collision or to end of chain */
        r[current_particle][0] += best_s * ex;
        r[current_particle][1] += best_s * ey;
        r[current_particle][2] += best_s * ez;

        wrap_coordinate(&r[current_particle][0], box[0]);
        wrap_coordinate(&r[current_particle][1], box[1]);
        wrap_coordinate(&r[current_particle][2], box[2]);

        /* Keep the linked-cell structure correct after this displacement */
        int new_cell = get_cell_index(r[current_particle][0],
                                      r[current_particle][1],
                                      r[current_particle][2]);
        update_cell_list_local(current_particle, old_cell, new_cell);

        remaining -= best_s;

        /* No collision before chain ends */
        if(hit_particle < 0){
            break;
        }

        /* Continue chain with the particle that was hit */
        current_particle = hit_particle;
    }

    return 1;
}


void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coord_steps/coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;

    build_cell_list();
}


double current_volume(void){
    double V = 1.0;
    for(int d = 0; d < NDIM; d++){
        V *= box[d];
    }
    return V;
}

void calc_mean_var(const double *data, int start, int end,
                   double *mean, double *var){
    int count = end - start;

    if(count <= 0){
        *mean = 0.0;
        *var = 0.0;
        return;
    }

    double sum = 0.0;
    for(int i = start; i < end; i++){
        sum += data[i];
    }

    *mean = sum / count;

    if(count <= 1){
        *var = 0.0;
        return;
    }

    double sq = 0.0;
    for(int i = start; i < end; i++){
        double diff = data[i] - *mean;
        sq += diff * diff;
    }

    *var = sq / (count - 1);
}

/* 
   Returns 1 if the volume appears equilibrated.
   We compare the last two volume blocks:
     old block:    [n - 2W, n - W)
     recent block: [n - W, n)

   Equilibrium criterion:
     |mean_recent - mean_old| < max(eq_rel_tol * mean_recent, 2 * standard_error)

   The relative tolerance catches slow drift.
   The standard-error tolerance prevents false rejection caused by normal NPT volume noise.
*/
int volume_is_equilibrated(const double *volume_history, int n_samples,
                           double *old_mean_out, double *new_mean_out,
                           double *diff_out, double *tol_out){
    int W = eq_block_samples;

    if(n_samples < 2 * W){
        return 0;
    }

    int start_old = n_samples - 2 * W;
    int end_old   = n_samples - W;
    int start_new = n_samples - W;
    int end_new   = n_samples;

    double mean_old, var_old;
    double mean_new, var_new;

    calc_mean_var(volume_history, start_old, end_old, &mean_old, &var_old);
    calc_mean_var(volume_history, start_new, end_new, &mean_new, &var_new);

    double diff = fabs(mean_new - mean_old);
    double standard_error = sqrt(var_old / W + var_new / W);
    double tol = eq_rel_tol * fabs(mean_new);

    if(tol < 2.0 * standard_error){
        tol = 2.0 * standard_error;
    }

    if(old_mean_out) *old_mean_out = mean_old;
    if(new_mean_out) *new_mean_out = mean_new;
    if(diff_out) *diff_out = diff;
    if(tol_out) *tol_out = tol;

    return (diff < tol);
}



void do_one_mc_step(int *move_accepted, int *move_trials,
                    int *vol_accepted, int *vol_trials){
    if (get_random_double(0.0, 1.0) < 0.90) {
        /* 90% probability: one local sweep */
        for (int n = 0; n < n_particles; ++n) {
            *move_accepted += move_particle_classic();
            *move_trials += 1;
        }
    }
    else {
        /* 10% probability: one event-chain move */
        *move_accepted += move_particle_chain();
        *move_trials += 1;
    }

    *vol_accepted += change_volume();
    *vol_trials += 1;
}

void append_volume_sample(double **history, int *n_samples, int *capacity, double V){
    if(*n_samples >= *capacity){
        int new_capacity = (*capacity == 0) ? 1024 : 2 * (*capacity);
        double *new_history = (double*)realloc(*history, new_capacity * sizeof(double));

        if(new_history == NULL){
            printf("Error: could not grow volume history array.\n");
            free(*history);
            exit(1);
        }

        *history = new_history;
        *capacity = new_capacity;
    }

    (*history)[*n_samples] = V;
    *n_samples += 1;
}

int main(int argc, char* argv[]){

    //timer
    clock_t start_time = clock();
    clock_t end_time;
    double cpu_time_used;

    time_t start_wall = time(NULL);

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();
    allocate_cell_list();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_packing_fraction();
    dsfmt_seed(time(NULL));

    /* Output file with one row per pressure. The reported volume is the
       mean of the last equilibrium block, not the mean over a fixed run. */
    FILE *csv_file = fopen("pVdata_until_equilibrium.csv", "w");
    if (csv_file == NULL) {
        printf("Error opening CSV file for writing.\n");
        return 1;
    }

    /* Trace file for plotting/checking the equilibration decision. */
    FILE *trace_file = fopen("equilibration_trace_until_equilibrium.csv", "w");
    if (trace_file == NULL) {
        printf("Error opening equilibration trace file for writing.\n");
        fclose(csv_file);
        return 1;
    }

    fprintf(csv_file, "Simulated NPT Data - run until equilibrium\n");
    fprintf(csv_file, "P_Solid,V_Solid,P_Liquid,V_Liquid,EquilibrationStep,EquilibriumDetected,BlockSamples\n");

    fprintf(trace_file, "betaP,step,volume,equilibrated,old_block_mean,new_block_mean,diff,tolerance,stable_checks\n");

    // Define the pressures for the liquid branch to check against Carnahan-Starling
    double liquid_pressures[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
   // double liquid_pressures[] = {0.5, 2.0, 3.0, 5.0, 8.0};

  //  double liquid_pressures[] = {3.0, 4.0, 5.0, 6.0};

    int num_pressures = sizeof(liquid_pressures) / sizeof(liquid_pressures[0]);

    /* Outer loop for pressure sweep.
       For every pressure, the code keeps running until the volume has stopped
       drifting according to volume_is_equilibrated(). */
    for (int p_idx = 0; p_idx < num_pressures; p_idx++) {
        betaP = liquid_pressures[p_idx];
        printf("\n--- Starting run for betaP = %.2lf ---\n", betaP);
        printf("Goal: run until equilibrium is detected, not until a fixed mc_steps value.\n");
        printf("#Step \t Volume \t Move-acc \t Vol-acc \t Stable checks \n");

        int move_accepted = 0;
        int move_trials = 0;
        int vol_accepted = 0;
        int vol_trials = 0;

        int equilibrated = 0;
        int equilibration_step = -1;
        int stable_checks = 0;

        double equilibrium_volume = 0.0;
        double old_mean = 0.0;
        double new_mean = 0.0;
        double diff = 0.0;
        double tol = 0.0;

        double *volume_history = NULL;
        int volume_samples = 0;
        int volume_capacity = 0;

        int step = 0;

        while(!equilibrated && step < max_safety_steps){
            step++;

            do_one_mc_step(&move_accepted, &move_trials,
                           &vol_accepted, &vol_trials);

            if(step % output_steps == 0){
                double V = current_volume();
                append_volume_sample(&volume_history, &volume_samples, &volume_capacity, V);

                if(step >= eq_min_steps){
                    if(volume_is_equilibrated(volume_history, volume_samples,
                                              &old_mean, &new_mean, &diff, &tol)){
                        stable_checks++;
                    }
                    else{
                        stable_checks = 0;
                    }

                    if(stable_checks >= eq_required_stable){
                        equilibrated = 1;
                        equilibration_step = step;
                        equilibrium_volume = new_mean;

                        printf("Equilibrium detected at step %d for betaP=%.2lf. "
                               "Recent block mean volume = %.6lf, previous block mean = %.6lf, "
                               "diff = %.6lf, tolerance = %.6lf\n",
                               step, betaP, new_mean, old_mean, diff, tol);
                    }
                }

                fprintf(trace_file, "%lf,%d,%lf,%d,%lf,%lf,%lf,%lf,%d\n",
                        betaP, step, V, equilibrated,
                        old_mean, new_mean, diff, tol, stable_checks);
            }

            if(step % progress_print_steps == 0){
                double move_acc_rate = (move_trials > 0) ? (double)move_accepted / move_trials : 0.0;
                double vol_acc_rate  = (vol_trials > 0) ? (double)vol_accepted / vol_trials : 0.0;

                printf("%d \t %lf \t %lf \t %lf \t %d \n",
                       step, current_volume(), move_acc_rate, vol_acc_rate, stable_checks);

                move_accepted = 0;
                move_trials = 0;
                vol_accepted = 0;
                vol_trials = 0;
            }
        }

        if(!equilibrated){
            printf("Warning: equilibrium was not detected for betaP=%.2lf before the safety limit of %d steps.\n",
                   betaP, max_safety_steps);
            printf("No trusted equilibrium volume is reported for this pressure. Try increasing max_safety_steps, "
                   "loosening eq_rel_tol slightly, or checking the volume trace.\n");

            if(volume_samples >= eq_block_samples){
                int start = volume_samples - eq_block_samples;
                double var_dummy;
                calc_mean_var(volume_history, start, volume_samples, &equilibrium_volume, &var_dummy);
            }
            else{
                equilibrium_volume = current_volume();
            }
        }

        printf("Result for betaP=%.2lf: V_equilibrium = %lf, equilibration_step = %d, detected = %d\n",
               betaP, equilibrium_volume, equilibration_step, equilibrated);

        fprintf(csv_file, "0.0,0.0,%lf,%lf,%d,%d,%d\n",
                betaP, equilibrium_volume, equilibration_step, equilibrated, eq_block_samples);

        free(volume_history);
    }

    fclose(csv_file);
    fclose(trace_file);

    printf("\nData successfully written to pVdata_until_equilibrium.csv!\n");
    printf("Equilibration trace written to equilibration_trace_until_equilibrium.csv!\n");

    //End bit of timing code
    end_time = clock();
    time_t end_wall = time(NULL);

    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nTotal (CPU) simulation time: %f seconds\n", cpu_time_used);
    printf("Wall clock time: %ld seconds\n", end_wall - start_wall);

    free(head);
    free(lscl);

    return 0;
}
