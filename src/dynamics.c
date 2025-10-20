#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"


// Calculate the total force on each particle
void forces_from_pos(Par *par, double *pos, double *force){
  int i, j, d;
  double r2;			// Distance squared between a pair of particles
  double f[D], dist[D];		// Force and distance between a pair of particles
  double *ipos, *jpos;		// Pointers to positions of particles i and j
  
  // Initialize the "total force" array to zero
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      force[D * i + d] = 0.0;

  // Loop over all pairs of particles
  for (i = 0; i < par->n - 1; i++) {
    ipos = pos + D * i;
    for (j = i+1; j < par->n; j++) {
      jpos = pos + D * j;
      
      r2 = dist2(par->L, ipos, jpos, dist);
      if (r2 < CUT * CUT) {
				// Each pair of interacting particles should come here
				// Calculate the force due do this interaction
				one_force(f, r2, dist);		// Calculate the vector f on the basis of r2 and dist
				for (d = 0; d < D; d++) {
					force[D * i + d] += f[d];
					force[D * j + d] -= f[d];
				}
      }
    }
  }
}

#ifdef BROWN
void pos_from_force(Par *par, double *pos, double *force){
	double c_brown = sqrt(6*par->T/(par->alpha*par->deltat));

	for (int i = 0; i < par->n; i++){
		for (int d = 0; d < D; d++){
			double xi = dran_sign();
			double newpos = pos[D * i + d] + (force[D*i+d] * par->deltat)/par->alpha + xi*c_brown*par->deltat;
			pos[D * i + d] = inrange(newpos, par->L[d]);
		}
	}
}
#endif

#ifdef VEL
// This function takes care of the Langevin terms by adding a term to the force:
// -alpha * velocity + noise
void langevin_forces(Par *par, double *vel, double *force)
{
  int i, d;
  double clang = sqrt(6*par->alpha*par->T/par->deltat);

  for (i = 0; i < par->n; i++){
    for (d = 0; d < D; d++) {
      // Fix this (5). Use dran_sign() which returns a value between -1 and 1.
      // Fix this (5):   force[D * i + d] += ???;
			double xi = dran_sign();
			force[D * i + d] += -par->alpha*vel[D * i + d] + xi*clang;
  	}
	}
}


// Here we use Newton's second law. For each particle, i, and direction, d, use
// vel(t + deltat) = vel(t) + (force / mass) * deltat
void vel_from_force(Par *par, double *vel, double *force)
{
  int i, d;

  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
   		vel[D * i + d] += 1/par->mass * force[D * i + d] * par->deltat;
    }
  }
}


// Calculate new positions for all particles. For each particle and direction, use
// pos(t + deltat) = pos(t) + vel * deltat
// Be careful to stay within the simulation box.
void pos_from_vel(Par *par, double *pos, double *vel)
{
  int i, d;
  
  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
      double newpos = pos[D * i + d] + par->deltat * vel[D * i + d];
      pos[D * i + d] = inrange(newpos, par->L[d]);
    }
  }
}
#endif 


// The function 'step' propagates the equations the time deltat.
// Three different kinds of dynamics:
// (1) Molecular dynamics (this is what we get when alpha = 0).
// (2) Langevin dynamics (noise and a damping term are added in 'langevin_forces')
// (3) Brownian dynamics (with no velocity term)
int step(Par *par, double *atoms, double *force)
{
  double *pos = atoms;

  // First determine the forces from the positions of all particles
  forces_from_pos(par, pos, force);

#ifdef VEL
  double *vel = atoms + D * par->n;
  // If Langevin simulations then add noise and damping to the forces
  if (par->alpha > 0.0)
    langevin_forces(par, vel, force);

  vel_from_force(par, vel, force);
  
  pos_from_vel(par, pos, vel);
#endif

#ifdef BROWN
  // Brownian dynamics goes here
  pos_from_force(par, pos, force);
#endif
  
  return 0;
}
