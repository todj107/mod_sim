#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"
extern int nnsize;


int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double r2, ediff;
  double newpos[D], dist[D];

#ifdef FAST
  check_neighbor_list(par, pos);
#endif
  
  for (i = 0; i < par->n; i++) {
    // Suggested new position for pos[i]
    for (d = 0; d < D; d++)
      newpos[d] = inrange(pos[D * i + d] + par->b * dran_sign(), par->L[d]);

    // Calculate the change in energy by moving particle i from pos[i] to newpos.
    ediff = 0.0;
// #ifdef FAST
//     int jj;
//     for (jj = 0; jj < neighbor_count[i]; jj++) {
//       j = neighbor_list[nnsize * i + jj];
// #else
		for (j = 0; j < par->n; j++) {
			if (j == i) continue;
// #endif
      
			// Note: "pos + D * i" is the same as "&pos[D * i]"
			r2 = dist2(par->L, pos + D * i, pos + D * j, dist);	// Distance^2 between particle i and j
			ediff -= pair_interaction(r2);

			r2 = dist2(par->L, newpos, pos + D * j, dist);	// Distance^2 between the new pos and part j
			ediff += pair_interaction(r2);
		}

		// Accept this new position?
		if ((ediff < 0.0) || exp(-ediff / par->T) > dran()) {
			for (d = 0; d < D; d++)
				pos[D * i + d] = newpos[d];
			naccept++;
		}
	}
	return naccept;
}
