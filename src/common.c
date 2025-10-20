#include "define.h"
#include "sim.h"
#include <math.h>

// This file includes some basic functions for calculating particle distances and forces.

#define ICUT6 (1.0 / (CUT * CUT * CUT * CUT * CUT * CUT))

// Function 'inrange' helps to implement periodic boundary conditions.
// Input: r, which has to be in the interval [-max, 2 max)
// Return value: r, in the interval [0,max).
double inrange(double r, double max)
{    
  if (r < 0.0)
    r += max;
  else
    if (r > max)
      r -= max;
  return r;
}


// *** pair_interaction ***
// Lennard-Jones interaction from Eq (1) in the lab instructions.
// Input: distance squared between two particles
// Return value: interaction energy
double pair_interaction(double r2)
{
  if (r2 > CUT * CUT)
    return 0.0;
	
	double u_r = 4*(1/pow(r2,6) - 1/pow(r2,3));
	double u_rc = 4*(ICUT6*ICUT6 - ICUT6);

	return u_r - u_rc;
}

// The magnitude of the force from Eq (2).
double force_magnitude(double r2)
{
	double r = sqrt(r2);
	double force;
	force = 1/r*(48/pow(r,12) - 24/pow(r,6));
	return force;
}


// Force determined with Eq (3).
// Input:  r2 = distance squared between two particles.
//         dist = vector with distances (x, y, ...)
// Output: f = vector with the force
// Note that this should work for general dimensions.
// We use 'D' for the dimensionality.
void one_force(double *f, double r2, double *dist)
{
  int d;
  double force = force_magnitude(r2);
  // Fix this (2): split the force up in the different directions
	for (int d = 0; d < D; d++){
		f[d] = force*dist[d]/sqrt(r2);
	}
}


// One dimensional distance (to the closest mirror point)
double distance(double L, double r1, double r2)
{
  double dist = r1 - r2;
  if (fabs(dist) > L / 2) {
    if (dist > 0.0)
      dist -= L;
    else 
      dist += L;
  }
  return dist;
}


// The function returns the D-dimensional distance squared.
// Input: L, p1, p2, which are all D-dimensional vectors.
//   L is the size of the simulation cell,
//   p1 and p2 are the particle positions,
// Output: the D-dimensional vector dist
// Return value: distance squared
double dist2(double *L, double *p1, double *p2, double *dist)
{
  int d;
  double sumdist2 = 0.0;

  for (d = 0; d < D; d++) {
    dist[d] = distance(L[d], p1[d], p2[d]);
    sumdist2 += dist[d] * dist[d];
  }

  return sumdist2;
}


