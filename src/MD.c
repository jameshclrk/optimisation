/*
*  Simple molecular dynamics code.
*  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
*
* This program implements:
*     long range inverse square forces between particles. F = G * m1*m2 / r**2
*     viscosity term     F = -u V
* If 2 particles approach closer than Size we flip the direction of the
* interaction force to approximate a collision.
*
* Coordinates are relative to a large central mass and the entire system is moving relative to the
* viscous media.
* If 2 particles approach closer than Size we flip the direction of the
* interaction force to approximate a collision.
*
* This program was developed as part of a code optimisation course
* and is therefore deliberately inefficient.
*
*/
#include <stdio.h>
#include <math.h>
#include "coord.h"

/* Helper macros to calculate the force */
#define TMP_FORCE(m1, m2, r, idx) G*m1*m2/(pow(r[idx],3));
#define FORCE_DIM(pos, idx, dim) tmp_force*pos[idx][dim]

/* Helper macro to calculate vector separation */
#define VECTOR_LENGTH(r, pos, i) r[i] = 0.0; \
for(d = 0; d < Ndim; d++) r[i] += (pos[i][d] * pos[i][d]); \
r[i] = sqrt(r[i]);

void evolve(int count,double dt){
    int step;
    int i, j, k, d;
    int collided;
    double m_size, tmp_force, tmp_force_2, dist;
    /*
    * Loop over timesteps.
    */
    for(step = 1; step <= count; step++){
        printf("timestep %d\n",step);
        printf("collisions %d\n",collisions);

        for(i = 0; i < Nbody; i++){
            /* calculate distance from central mass */
            VECTOR_LENGTH(r, pos, i)

            /* calculate central force */
            /* set the viscosity term in the force calculation */
            /* add the wind term in the force calculation */
            tmp_force = TMP_FORCE(mass[i], M_central, r, i);
            for(d = 0; d < Ndim; d++) {
                f[i][d] = -visc[i] * (vel[i][d] + wind[d]) - FORCE_DIM(pos, i, d);
            }
        }
        /* calculate pairwise separation of particles */
        for(i = 0, k = 0; i < Nbody; i++){
            for(j = i+1; j < Nbody; j++, k++){
                for(d = 0; d < Ndim; d++){
                    delta_pos[k][d] = pos[i][d] - pos[j][d];
                }
            }
        }
        /* calculate norm of seperation vector */
        for(i = 0; i < Npair; i++){
            VECTOR_LENGTH(delta_r, delta_pos, i)
        }
        /*
        * add pairwise forces.
        */
        for(i = 0, k = 0; i < Nbody; i++){
            for(j = i + 1; j < Nbody; j++, k++){
                m_size = radius[i] + radius[j];
                /* Calculate a temporay force variable */
                tmp_force = TMP_FORCE(mass[i], mass[j], delta_r, k);
                /*  flip force if close in */
                if (delta_r[k] < m_size ) {
                    tmp_force *=-1.0;
                    collisions++;
               }
                for(d = 0; d < Ndim; d++){
                    /* Use the temporary force variable to calculate the real force */
                    tmp_force_2 = FORCE_DIM(delta_pos, k, d);
                    f[i][d] = f[i][d] - tmp_force_2;
                    f[j][d] = f[j][d] + tmp_force_2;
                }
            }
        }

        for(i = 0; i < Nbody; i++){
            for(d = 0; d < Ndim; d++){
               /* update positions */
                pos[i][d] = pos[i][d] + dt * vel[i][d];
                /* update velocities */
                vel[i][d] = vel[i][d] + dt * (f[i][d]/mass[i]);
            }
        }
    }
}
