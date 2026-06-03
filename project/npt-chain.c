#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 500

/* Initialization variables */
const int mc_steps = 2000;
const int output_steps = 10;
const double packing_fraction = 0.3;
const double diameter = 1.0;
const double delta  = 0.1;
/* Volume change -deltaV, delta V */
const double deltaV = 2.0;
/* Reduced pressure \beta P */
const char* init_filename = "xyz_500.dat"; //<----------------------aanpassen

/* Simulation variables */
double betaP;
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];
double chain_len = 10.; //<--------------------------- nieuwe var

int change_volume(void){

    double Vold = box[0]*box[1]*box[2];
    double dV = deltaV * (2.0*dsfmt_genrand() - 1.0);
    double Vnew = Vold + dV;
    if (Vnew <= 0.0) return 0;

    double s = cbrt(Vnew / Vold);
    double box_old[NDIM];


    for(int d = 0; d < NDIM; d++) box_old[d] = box[d];

    static double r_old[N][NDIM];
    for(int i = 0; i < n_particles; i++){
        for(int d = 0; d < NDIM; d++){
            r_old[i][d] = r[i][d];
        }
    }

    for(int i = 0; i < n_particles; i++){
        for(int d = 0; d < NDIM; d++){
            r[i][d] *= s;
        }
    }
    for(int d = 0; d < NDIM; d++){
        box[d] *= s;
    }

    for(int i = 0; i < n_particles; i++){
        for(int j = i + 1; j < n_particles; j++){
            double dx = r[i][0] - r[j][0];
            double dy = r[i][1] - r[j][1];
            double dz = r[i][2] - r[j][2];

            if (dx >  0.5*box[0]) dx -= box[0];
            if (dx < -0.5*box[0]) dx += box[0];
            if (dy >  0.5*box[1]) dy -= box[1];
            if (dy < -0.5*box[1]) dy += box[1];
            if (dz >  0.5*box[2]) dz -= box[2];
            if (dz < -0.5*box[2]) dz += box[2];

            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < diameter*diameter){
                /* reject: revert everything */
                for(int d = 0; d < NDIM; d++) box[d] = box_old[d];
                for(int a = 0; a < n_particles; a++){
                    for(int d = 0; d < NDIM; d++){
                        r[a][d] = r_old[a][d];
                    }
                }
                return 0;
            }
        }
    }

    double exponent = -betaP * (Vnew - Vold) + (double)n_particles * log(Vnew / Vold);

    double u = dsfmt_genrand();
    if (log(u) < exponent){return 1;}
    else{
        for(int d = 0; d < NDIM; d++) {
            box[d] = box_old[d];
        }
        for(int a = 0; a < n_particles; a++){
            for(int d = 0; d < NDIM; d++){
                r[a][d] = r_old[a][d];
            }
        }
        return 0;
    }
}

void read_data(void){
    FILE *fptr;
    fptr = fopen(init_filename, "r");


    //number of particles
    int num;
    fscanf(fptr, "%d", &num);
    n_particles = num;
    printf("n particles: %d\n", n_particles);


    //checking len for all coords
    float x1, x2, y1, y2, z1, z2;
    fscanf(fptr, "%f %f", &x1, &x2);
    fscanf(fptr, "%f %f", &y1, &y2);
    fscanf(fptr, "%f %f", &z1, &z2);

    box[0] = (double)x2;
    box[1] = (double)y2;
    box[2] = (double)z2;


    //coords into array and making a fourth var d to read it but not use
    float x, y, z, d;
    for (int i = 0; i < N; i++) {
        fscanf(fptr, "%f %f %f %f", &x, &y, &z, &d);
        r[i][0] = (double)x;
        r[i][1] = (double)y;
        r[i][2] = (double)z;
    }

    radius = d/2; 
    fclose(fptr);
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

    /* Pick random starting particle.
       Your old code used n_particles - 1, which never selected the last particle. */
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

        /* Find nearest collision in the chain direction */
        for(int j = 0; j < n_particles; j++){

            if(j == current_particle) continue;

            double dx = x - r[j][0];
            double dy = y - r[j][1];
            double dz = z - r[j][2];

            dx = min_image(dx, box[0]);
            dy = min_image(dy, box[1]);
            dz = min_image(dz, box[2]);

            /*
               Collision condition:
               |rij + s e|^2 = sigma^2

               b = rij . e
               c = rij^2 - sigma^2
               s = -b - sqrt(b^2 - c)
            */

            double b = dx*ex + dy*ey + dz*ez;

            /* If b >= 0, particle j is not in front of the moving particle */
            if(b >= 0.0) continue;

            double c = dx*dx + dy*dy + dz*dz - sigma2;
            double discriminant = b*b - c;

            if(discriminant <= 0.0) continue;

            double s = -b - sqrt(discriminant);

            if(s > eps && s < best_s){
                best_s = s;
                hit_particle = j;
            }
        }

        /* Move current particle either to next collision or to end of chain */
        r[current_particle][0] += best_s * ex;
        r[current_particle][1] += best_s * ey;
        r[current_particle][2] += best_s * ez;

        wrap_coordinate(&r[current_particle][0], box[0]);
        wrap_coordinate(&r[current_particle][1], box[1]);
        wrap_coordinate(&r[current_particle][2], box[2]);

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


/*
int move_particle_chain(void){

    double picker = (double)(n_particles - 1)*dsfmt_genrand();
    double chain_len_local = chain_len;
    int randparticle = (int)picker; 
    int current_particle = randparticle;
    int new_particle;

    double x = r[randparticle][0];
    double y = r[randparticle][1];
    double z = r[randparticle][2];

    double deltax = delta*(2*dsfmt_genrand() - 1);
    double deltay = delta*(2*dsfmt_genrand() - 1);
    double deltaz = delta*(2*dsfmt_genrand() - 1);
    double deltar = sqrt(pow(deltax,2) + pow(deltay,2) + pow(deltaz,2));

    //in bepaalde richting worden alle particles met een chain reaction gestuurd
    while(chain_len_local > 0){
        new_particle = move_particle(current_particle, x, y, z, deltax, deltay, deltaz);
        current_particle = new_particle;
        chain_len_local -= deltar;
        return 1;
    }

    return 0;
}*/

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.dat", step);
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
    //allocate_cell_list();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_packing_fraction();
    dsfmt_seed(time(NULL));

    // Open CSV for Python script
    FILE *csv_file = fopen("pVdata.csv", "w");
    if (csv_file == NULL) {
        printf("Error opening CSV file for writing.\n");
        return 1;
    }
    
    // Write two header lines because Python script uses data[2:,:]
    fprintf(csv_file, "Simulated NPT Data\n");
    fprintf(csv_file, "P_Solid,V_Solid,P_Liquid,V_Liquid\n");

    // Define the pressures for the liquid branch to check against Carnahan-Starling
    double liquid_pressures[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
   // double liquid_pressures[] = {3.0, 4.0, 5.0, 6.0};

    int num_pressures = sizeof(liquid_pressures) / sizeof(liquid_pressures[0]);

    // Outer loop for pressure sweep (Continuous compression)
    for (int p_idx = 0; p_idx < num_pressures; p_idx++) {
        betaP = liquid_pressures[p_idx];
        printf("\n--- Starting run for betaP = %.2lf ---\n", betaP);
        printf("#Step \t Volume \t Move-acc \t Vol-acc \n");

        int move_accepted = 0;
        int vol_accepted = 0;
        double mean_sum = 0.0;
        int samples = 0;

        for(int step = 0; step < mc_steps; ++step){
            for(int n = 0; n < n_particles; ++n){
                move_accepted += move_particle_chain();
            }
            vol_accepted += change_volume();

            // Safely compute mean volume without hardcoded arrays
            if(step > 100 && step % output_steps == 0){
                mean_sum += (box[0] * box[1] * box[2]);
                samples++;
            }

            if(step % (mc_steps / 10) == 0){ // Print to terminal 10 times per run
                printf("%d \t %lf \t %lf \t %lf \n", 
                        step, box[0] * box[1] * box[2], 
                        (double)move_accepted / (n_particles * (mc_steps/10)), 
                        (double)vol_accepted / (mc_steps/10));
                move_accepted = 0;
                vol_accepted = 0;
                // write_data(step); // Optional: Uncomment if you want to save coordinates
            }
        }

        double mean_volume = mean_sum / samples;
        printf("Result: Mean volume for betaP=%.2lf is %lf\n", betaP, mean_volume);

        // Write directly to CSV. Columns 0 and 1 are 0.0 since we are sweeping the liquid phase
        fprintf(csv_file, "0.0,0.0,%lf,%lf\n", betaP, mean_volume);
    }

    fclose(csv_file);
    printf("\nData successfully written to pVdata.csv!\n");

    return 0;
}
