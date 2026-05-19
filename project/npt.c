#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 500

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;
const double packing_fraction = 0.6;
const double diameter = 1.0;
const double delta  = 0.1;
/* Volume change -deltaV, delta V */
const double deltaV = 2.0;
/* Reduced pressure \beta P */
const double betaP = 15.0;
const char* init_filename = "xyz2.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];
double volume_list[900] = {}; //length is hard coded, because we have conditions on which volumes we take 
// We check volume every (output_steps) steps, but dont look at the first 10000 steps, when the systen is still equilibrating

//Self made variables for easier computing
double diameter_squared  = diameter*diameter;

//Generates random double between any max and min value
double get_random_double(double min, double max) {
    // Swap if min is accidentally greater than max
    if (min > max) {
        double temp = min;
        min = max;
        max = temp;
    }
    
    // Scale the 0.0 to 1.0 result to your desired range
    return min + (dsfmt_genrand() * (max - min));
}

//Applies periodic boundary conditions on a coordinate
double apply_pbc(double position, double box_length) {
    if (position >= box_length) position -= box_length;
    if (position < 0.0) position += box_length;
    return position;
}


/* Functions */
int change_volume(void){
    double dV = get_random_double(-1*deltaV, deltaV);
    
    //Calculate volume_old 
    double volume_old = 1.0;
    double box_old[NDIM];
    for (int i = 0; i < NDIM; i++){
        volume_old *= box[i];
        box_old[i] = box[i];
    }
    
    //Calculate new volume
    double volume_new = volume_old + dV;
    if (volume_new <= 0){
        return 0;
    }

    //Set new box dimensions
    double half_length[NDIM];
    for (int i = 0; i < NDIM; i++){
        box[i] = cbrt(volume_new);
        half_length[i] = 0.5 * box[i];
    }
    
    //Scale ALL particles to the new box
    double scale_x = box[0] / box_old[0];
    double scale_y = box[1] / box_old[1];
    double scale_z = box[2] / box_old[2];
    
    for (int n = 0; n < n_particles; n++){
        r[n][0] *= scale_x;
        r[n][1] *= scale_y;
        r[n][2] *= scale_z;
    }
    
    if (volume_new < volume_old){ //only need to check for overlap if volume is shrinking
        //Check for overlaps
        for (int n = 0; n < n_particles; n++){
            double rx = r[n][0];
            double ry = r[n][1];
            double rz = r[n][2];
            for (int m = n + 1; m < n_particles; m++){
                double dx = rx - r[m][0];
                double dy = ry - r[m][1];
                double dz = rz - r[m][2];
                
                //Apply nearest image convention for X
                if (dx > half_length[0]) dx -= box[0];
                if (dx <= -half_length[0]) dx += box[0];
                
                double dx2 = dx * dx; //shortcut, if dx^2 already bigger than diameter^2, no overlap possible
                if (dx2 >= diameter_squared) continue;

                //Apply nearest image convention for Y
                if (dy > half_length[1]) dy -= box[1];
                if (dy <= -half_length[1]) dy += box[1];
                
                double dy2 = dy * dy; //same shortcut now for x and y
                if (dx2 + dy2 >= diameter_squared) continue;

                //Apply nearest image convention for Z
                if (dz > half_length[2]) dz -= box[2];
                if (dz <= -half_length[2]) dz += box[2];
                
                double dz2 = dz * dz;
                
                //Final check uses the cached values
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
    

    //calculate acceptance probability
    double acc_arg = -(betaP / (diameter_squared * diameter)) * (volume_new - volume_old) + n_particles * log(volume_new / volume_old);
    double acc_exp = exp(acc_arg);
    double acceptance = fmin(1.0, acc_exp);
    
    //Accept or Reject
    if(dsfmt_genrand() < acceptance){
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
    //Opening file
    FILE *file = fopen(init_filename, "r");
    if (file == NULL){
        printf("File not found\n");
        exit(1);
    }

    //Reads number of particles
    if(fscanf(file, "%d", &n_particles) != 1){
        printf("%d\n",n_particles);
        fclose(file);
        exit(1);
    }


    //Test to see if n_particles is read correctly, uncomment to check
    //printf("Read n_particles: %d\n", n_particles);


    //Read dimensions of the box. Data file has coords of edges, but this saves it as lengths of edges
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

    //Reads positions of particles. If NDIM=2, forgets the z-component.
    double temp_diameter;
    if (NDIM == 2){
        for(int n = 0; n < N; n++){
        if (fscanf(file, "%lf %lf %*lf %lf", &r[n][0], &r[n][1], &temp_diameter) != 3) {
            printf("Error: Failed to read coordinates for particle %d.\n", n);
            fclose(file);
            exit(1);
            }
            //Tests if particle positions are shows. Uncomment to check
            //if (n < 5 || n == n_particles - 1) {
            //    printf("Particle %d: X=%lf, Y=%lf (Z skipped, Diameter=%lf)\n", n, r[n][0], r[n][1], temp_diameter);
            //}
        }
    }
    else if (NDIM == 3){
        for(int n = 0; n < N; n++){
        if (fscanf(file, "%lf %lf %lf %lf", &r[n][0], &r[n][1], &r[n][2], &temp_diameter) != 4) {
            printf("Error: Failed to read coordinates for particle %d.\n", n);
            fclose(file);
            exit(1);
            }
            //Tests if particle positions are shows. Uncomment to check
            //if (n < 5 || n == n_particles - 1) {
            //    printf("Particle %d: X=%lf, Y=%lf, Z=%lf (Diameter=%lf)\n", n, r[n][0], r[n][1], r[n][2], temp_diameter);
            //}
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

int move_particle(void){
    //Get particle index and displacement
    int particle_index = (int) get_random_double(0,n_particles);
    double half_length[NDIM];
    for (int i = 0; i < NDIM; i++){
        half_length[i] = 0.5*box[i];
    }

    //save old r incase displacement is rejected
    double old_r[NDIM];
    for (int i = 0; i < NDIM; i++){
        old_r[i] = r[particle_index][i];
    }

    //For each index add displacement and apply periodic boundary conditions
    for (int i = 0; i < NDIM; i++){
        r[particle_index][i] = r[particle_index][i] + get_random_double(-1*delta, delta);
        r[particle_index][i] = apply_pbc(r[particle_index][i], box[i]); //this function is defined above
    }
    double rx = r[particle_index][0];
    double ry = r[particle_index][1];
    double rz = r[particle_index][2];
    for (int n = 0; n < n_particles; n++){
            if (particle_index != n){
                double dx = rx - r[n][0];
                double dy = ry - r[n][1];
                double dz = rz - r[n][2];
                //Apply nearest image convention for all directions
                if (dx > half_length[0]){
                    dx = dx - box[0];
                }
                 if (dx <= -half_length[0]){
                    dx = dx + box[0];
                }
                //If x distance is already greater that dist^2, it cant overlap, so this is a loop shortcut to save time
                if (dx * dx >= diameter_squared) continue;

                if (dy > half_length[1]){
                    dy = dy - box[1];
                }
                 if (dy <= -half_length[1]){
                    dy = dy + box[1];
                }
                //Same shortcut for y and x together
                if (dy * dy + dx * dx >= diameter_squared) continue;

                if (dz > half_length[2]){
                    dz = dz - box[2];
                }
                 if (dz <= -half_length[2]){
                    dz = dz + box[2];
                }
                //Check if particles overlap, if they do, set r back to old values
                double distance_squared = dx*dx + dy*dy + dz*dz;
                if (distance_squared < diameter_squared){
                    for(int i = 0; i < NDIM; i++){
                        r[particle_index][i] = old_r[i];
                    }
                    return 0;
                }
            }
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
}

int main(int argc, char* argv[]){

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_packing_fraction();

    dsfmt_seed(time(NULL));
            
    printf("#Step \t Volume \t Move-acceptance\t Volume-acceptance \n");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            move_accepted += move_particle();
        }
        vol_accepted += change_volume();

        if(step%output_steps == 0 && step > 10000){
        volume_list[step/output_steps] = box[0]*box[1]*box[2];
        }

        if(step % output_steps == 0){
            printf("%d \t %lf \t %lf \t %lf \n", 
                    step, box[0] * box[1] * box[2], 
                    (double)move_accepted / (n_particles * output_steps), 
                    (double)vol_accepted /  output_steps);
            move_accepted = 0;
            vol_accepted = 0;
            write_data(step);
        }
    }

    double mean_sum;
    double mean_volume;
    for (int i = 0; i < mc_steps/output_steps; i++){
        mean_sum += volume_list[i];
    }
    mean_volume = mean_sum/(900); //hardcoded based on volume_list length
    printf("The mean volume is %lf\n", mean_volume);

    return 0;
}
