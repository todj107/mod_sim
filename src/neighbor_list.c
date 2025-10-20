#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
#define NNSIZE 80

#ifdef FAST
int *neighbor_list = NULL, *neighbor_count;
double *save_pos;
int nnsize = NNSIZE;


void update_neighbor_list(Par *par, double *pos)
{
  int i, j, id, count;
  double dist[D];
  if (!neighbor_list) {
    neighbor_count = malloc(par->n * sizeof(int));
    neighbor_list = malloc(nnsize * par->n * sizeof(int));
    save_pos = malloc(D * par->n * sizeof(double));
  }

  ...
}

void check_neighbor_list(Par *par, double *pos)
{
  ...
}
#endif
