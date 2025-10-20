#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran.h"

#define FNAMESIZE 100

typedef struct Measured {
  double epot, etot, pressure, etot2;
  // Fix this (4). Introduce variable for total energy squared
#ifdef VEL
  double ekin;
#endif
} Measured;

// All the values that define the model and the properties
// of the run are stored in a struct variable of type "Par".
typedef struct Par {
  int nblock, nsamp, ntherm, seed, n, df;
  double L[D], vol, rho, mass;
  double T;
#ifdef MC
  double b;
#else
  double alpha;
  double deltat;
#endif
} Par;

// In measure.c:
extern void m_initialize(Par *par);
extern void measure(Par *par, double *atoms, Measured *val);
extern void m_energy_print(FILE *estream, int is);
extern void m_block_averages(Par *par);
extern void m_print_results(Par *par);
extern void m_print_recent_values(Par *par);

// In sim.c:
extern char *add_strings(char *str1, char *str2);
extern char *get_filename(Par *par);

// In common.c:
extern double pair_interaction(double r2);
extern double force_magnitude(double r2);
extern void one_force(double *f, double r2, double *dist);
extern double inrange(double r, double max);
extern double distance(double L, double r1, double r2);
extern double dist2(double *L, double *p1, double *p2, double *dist);
extern void forces_from_pos(Par *par, double *pos, double *force);
#ifdef BROWN 
extern void pos_from_force(Par *par, double *pos, double *force);
#endif
#ifdef VEL
extern void langevin_forces(Par *par, double *vel, double *force);
extern void vel_from_forces(Par *par, double *vel, double *force);
extern void pos_from_vel(Par *par, double *pos, double *vel);
#endif
extern int step(Par *par, double *atoms, double *force);

// In config.c:
extern void init_vel(Par *par, double *vel);
extern void init_pos(Par *par, double *pos);
extern double *init_conf(Par *par);
extern double *read_conf(Par *par, double *atoms, char *fname, int just_attempting);
extern int write_conf(Par *par, double *atoms, char *fname);

#ifdef MC
// In (new file) monte_carlo.c?
extern int monte_carlo_sweep(Par *par, double *pos);
extern void check_neighbor_list(Par *par, double *pos);
#endif


#ifdef FAST
extern int *neighbor_list, *neighbor_count;
extern void check_neighbor_list(Par *par, double *pos);
#endif
  
