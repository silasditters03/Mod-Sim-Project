#include <stdio.h>
#include <math.h>
#include <time.h>

double apply_pbc(double position, double box_length) {
    // fmod gets the remainder of floating-point division
    double wrapped = fmod(position, box_length);
    
    // fmod can return a negative number if position is negative.
    // We add the box_length to wrap it correctly back into the positive domain.
    if (wrapped < 0.0) {
        wrapped += box_length;
    }
    
    return wrapped;
}

int main(){
    //Starting parameters
    int Nx = 4;
    int Ny = 4;
    int Nz = 4;
    int N = 4 * Nx * Ny * Nz; //Because 4 particles per cell
    double diameter = 1;
    double spacing = 2;

    //Dimensions box
    double Lx = Nx * spacing;
    double Ly = Ny * spacing;
    double Lz = Nz * spacing;

    // Define origin shifts (Set to 0.0 to place on the edges)
    // You can change these to shift the entire lattice by any arbitrary amount
    double shift_x = 0.0; 
    double shift_y = 0.0; 
    double shift_z = 0.0;

    // 2. Define the 4 basis vectors for an FCC unit cell
    // These are fractions of the lattice spacing 'l'
    double basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    };

    FILE *file = fopen("xyz2.dat", "w");

    // 3. Write the header
    fprintf(file, "%d\n", N);
    fprintf(file, "%f %f\n", 0, Lx);
    fprintf(file, "%f %f\n", 0, Ly);
    fprintf(file, "%f %f\n", 0, Lz);

    // 4. Generate the lattice and apply PBC
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {

                for(int b = 0; b < 4; b++){

                // Calculate raw positions with your arbitrary shift
                double x_raw = ((i + basis[b][0] )* spacing) + shift_x;
                double y_raw = ((j + basis[b][1] ) * spacing) + shift_y;
                double z_raw = ((k + basis[b][2] ) * spacing) + shift_z;

                // Apply Periodic Boundary Conditions to wrap them into the box
                double x_pbc = apply_pbc(x_raw, Lx);
                double y_pbc = apply_pbc(y_raw, Ly);
                double z_pbc = apply_pbc(z_raw, Lz);

                // Write the strictly bounded coordinates and the diameter
                fprintf(file, "%f %f %f %f\n", x_pbc, y_pbc, z_pbc, diameter);

                }
                
              
            }
        }
    }

    fclose(file);

    return 0;
}

/*
Answer to 4d:
The maximum packing fraction of a face-centered cubic lattice with spacing = diameter (touching spheres) is the ratio 
of the number of particles time the volume per particle, divide by the volume of the cell.
P = N_particle V_particle / V_cell = 4/3 pi (d/2)^3 / (sqrt(2)*d)^3 = pi/(3*sqrt(2)) is about 0.74...

So the ratio P_cubic/P_fcc = pi/6 / 0.74 = 0.707...
So FCC has a higher packing fraction
*/