#include <stdio.h>
#include <math.h>
#include <time.h>

double apply_pbc(double position, double box_length) {
    // fmod gets the remainder of floating-point division
    double wrapped = fmod(position, box_length);
    
    // In C, fmod can return a negative number if position is negative.
    // We add the box_length to wrap it correctly back into the positive domain.
    if (wrapped < 0.0) {
        wrapped += box_length;
    }
    
    return wrapped;
}

int main(){
    //Starting parameters
    int Nx = 10;
    int Ny = 10;
    int Nz = 10;
    int N = Nx * Ny * Nz;
    double diameter = 1;
    double spacing = 2;

    //Dimensions box
    double Lx = Nx * spacing;
    double Ly = Ny * spacing;
    double Lz = Nz * spacing;

    // 2. Define origin shifts (Set to 0.0 to place on the edges)
    // You can change these to shift the entire lattice by any arbitrary amount
    double shift_x = 0.0; 
    double shift_y = 0.0; 
    double shift_z = 0.0;

    FILE *file = fopen("xyz1.dat", "w");

    // 3. Write the header
    fprintf(file, "%d\n", N);
    fprintf(file, "%f %f\n", 0, Lx);
    fprintf(file, "%f %f\n", 0, Ly);
    fprintf(file, "%f %f\n", 0, Lz);

    // 4. Generate the lattice and apply PBC
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                
                // Calculate raw positions with your arbitrary shift
                double x_raw = (i * spacing) + shift_x;
                double y_raw = (j * spacing) + shift_y;
                double z_raw = (k * spacing) + shift_z;

                // Apply Periodic Boundary Conditions to wrap them into the box
                double x_pbc = apply_pbc(x_raw, Lx);
                double y_pbc = apply_pbc(y_raw, Ly);
                double z_pbc = apply_pbc(z_raw, Lz);

                // Write the strictly bounded coordinates and the diameter
                fprintf(file, "%f %f %f %f\n", x_pbc, y_pbc, z_pbc, diameter);
            }
        }
    }

    fclose(file);

    return 0;
}

/*
Answer to 4b:
The maximum packing fraction of a cubic lattice with spacing = diameter (touching spheres) is the ratio of the box volume to
the sphere volume P = V_sphere / V_cube = 4/3 pi (d/2)^3 / d^3 = pi/6
*/