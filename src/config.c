#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"


// Initialize the velocity with random numbers to get mv^2/2 = T/2 for each degree of freedom
#ifdef VEL
void init_vel(Par *par, double *vel)
{
  int i, d;
  double tot[D];

  // Initialize the random number generator (if not already done)
  init_ran(par->seed);

  for (d = 0; d < D; d++)
    tot[d] = 0.0;

  // First generate random values with <vel^2> = 1
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      tot[d] += vel[D * i + d] = sqrt(3.0) * dran_sign();

  // Avoid center of mass motion
  for (d = 0; d < D; d++)
    tot[d] /= par->n;

  double sumv2 = 0.0;
  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
      vel[D * i + d] -= tot[d];
      sumv2 += vel[D * i + d] * vel[D * i + d];
    }
  }

  // We want sumv2 = D * N * T
  double fact = sqrt((D * par->n * par->T) / sumv2);
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      vel[D * i + d] *= fact;
}
#endif


// Here we put the atoms in the simulation box in a square or cubic lattice
// just to avoid getting particles too close to each other which would give very big forces. 
void init_pos(Par *par, double *pos)
{
  int i, d;
  int m = ceil(exp( log(par->n) / D ));

#if D == 2
  for (i = 0; i < par->n; i++) {
    pos[D * i] = i % m;
    pos[D * i + 1] = i / m;
  }
#endif
#if D == 3
  for (i = 0; i < par->n; i++) {
    pos[D * i] = i % m;
    pos[D * i + 1] = (i / m) % m;
    pos[D * i + 2] = i / (m * m);
  }
#endif

  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++)
      pos[D * i + d] = (pos[D * i + d] + 0.5) * par->L[d] / m;
  } 
}


double *init_conf(Par *par)
{
  printf("Generate a simple starting configuration\n");
  double *atoms;
  atoms = malloc(par->n * par->df * sizeof(double));

  double *pos = atoms;
  init_pos(par, pos);

#ifdef VEL
  double *vel = atoms + D * par->n;
  init_vel(par, vel);
#endif

  return atoms;
}



double *read_conf(Par *par, double *atoms, char *name, int just_attempting)
{
  // If null pointer (i.e. no file name) determine a file name from the parameters
  if (!name)
    name = get_filename(par);

  char hash;
  FILE *stream;
  int n0, df0;
  double rho0, t0, L0;
  int i, d, num;
  
  stream = fopen(add_strings("conf/", name), "r");
  if (!stream) {
    if (!just_attempting)
      printf("*** Cannot open conf/%s\n", name);
    return NULL;
  }

  // Read first line to get the geometry of the file

  char *lineptr = NULL;
  size_t linelen;
  ssize_t bytes = getline(&lineptr, &linelen, stream);
  //  printf("Written %ld char to address %x\n", linelen, lineptr);
  
  num = sscanf(lineptr, "%c %d %d  %lf  %lf  %lf", &hash, &n0, &df0, &rho0, &t0, &L0);
  free(lineptr);
  //  printf("File with N=%d  df=%d\n", n0, df0);

  // Consistency checks
  if (num < 3) {
    fprintf(stderr, "*** Cannot read number of particles and number of degrees of freedom from conf/%s\n", name);
    perror(NULL);
    return 0;
  }

  if (n0 != par->n) {
    fprintf(stderr, "*** %d particles in conf/%s; %d particles expected.\n", par->n, name, n0);
    fclose(stderr);
    return NULL;
  }


  printf("Read from conf/%s with data for %d particles\n", name, par->n);

  double *aread = malloc(par->n * df0 * sizeof(double));

  num = 0;
  for (i = 0; i < par->n; i++) {
    double *ia = aread + df0 * i;
    int j;
    for (j = 0; j < df0; j++)
      num += fscanf(stream, "%lf ", ia + j);
  }
  fclose(stream);

  if (num < df0 * par->n) {
    fprintf(stderr, "Only %d values for particle %d (expected %d)\n", num, i, df0);
    return 0;
  }

  
  if (!atoms)
    atoms = malloc(par->n * par->df * sizeof(double));

  // Copy from (the temporary) aread to atoms
  for (i = 0; i < par->n; i++) {
    double *ia = aread + df0 * i;
    double *ipos = atoms + D * i;
    for (d = 0; d < D; d++)
      ipos[d] = ia[d];
#ifdef VEL
    if (df0 > D) {	// If there is velocity data in the file
      double *vel = atoms + D * par->n;
      double *ivel = vel + D * i;
      for (d = 0; d < D; d++)
	ivel[d] = ia[D + d];
    }
#endif
  }

  free(aread);

#ifdef VEL
  if (par->df > df0) {
    double *vel = atoms + D * par->n;
    init_vel(par, vel);
  }
#endif

  return atoms;
}


int write_conf(Par *par, double *atoms, char *name)
{
  FILE *stream;
  int i, d;
  double *pos = atoms;
  
  stream = fopen(add_strings("conf/", name), "w");
  if (!stream) {
    fprintf(stderr, "*** Cannot open conf/%s for write ", name);
    return 0;
  }

  // First print number of particles and number of degrees of freedom
  fprintf(stream, "# %d  %d  %g  %g  %g\n", par->n, par->df, par->rho, par->T, par->L[0]);

  // Then print position and velocity
  for (i = 0; i < par->n; i++) {
    double *ipos = pos + D * i;
    for (d = 0; d < D; d++)
      fprintf(stream, "%g  ", ipos[d]);
#ifdef VEL
    double *vel = atoms + D * par->n;
    double *ivel = vel + D * i;
    for (d = 0; d < D; d++)
      fprintf(stream, "%g  ", ivel[d]);
#endif
    fprintf(stream, "\n");
  }
  fclose(stream);

  printf("Configuration with %d particles written to conf/%s\n", par->n, name);
  return 1;
}
