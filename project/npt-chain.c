#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;
const double packing_fraction = 0.6;
const double diameter = 1.0;
const double delta  = 0.1;
/* Volume change -deltaV, delta V */
const double deltaV = 2.0;
/* Reduced pressure \beta P */
const double betaP = 3.0;
const char* init_filename = "xyz.dat"; //<----------------------aanpassen

/* Simulation variables */
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
}

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
            move_accepted += move_particle_chain(); //<------------------------------------aangepast
        }
        vol_accepted += change_volume();

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

    return 0;
}
