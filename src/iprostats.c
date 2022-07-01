#include "lise.h"

void createiprostats(int num_atoms, IPRO_STATS *restrict ipro_st) {
    ipro_st->dists = (int **)(calloc(num_atoms, sizeof(int *)));
    for (int i = 0; i < num_atoms; i++) {
        ipro_st->dists[i] = (int *)(calloc(num_atoms, sizeof(int)));
    }
}
int setiprostats(const IPRO *restrict ipro, IPRO_STATS *restrict ipro_st) {
    int tmp;
    int ret = 0;

    initstats(&ipro_st->X);
    initstats(&ipro_st->Y);
    initstats(&ipro_st->Z);

    for (int i = 0; i < ipro->numAtoms; i++) {
        setstats(ipro->atoms[i].x, &ipro_st->X);
        setstats(ipro->atoms[i].y, &ipro_st->Y);
        setstats(ipro->atoms[i].z, &ipro_st->Z);

        for (int j = i+1; j < ipro->numAtoms; j++) {
            tmp = sqrt(pow(ipro->atoms[i].x - ipro->atoms[j].x, 2) + pow(ipro->atoms[i].y - ipro->atoms[j].y, 2) + pow(ipro->atoms[i].z - ipro->atoms[j].z, 2));
            
            if (tmp < 0)
                return ret;
			if (tmp < 2)
				ipro_st->dists[i][j] = 0;
			else if (tmp > 13)
				ipro_st->dists[i][j] = 0;
			else if (tmp > 10)
				ipro_st->dists[i][j] = 2;
			else
				ipro_st->dists[i][j] = 1;
            
            ret++;
        }
    }

    finstats(&ipro_st->X);
    finstats(&ipro_st->Y);
    finstats(&ipro_st->Z);

    return ret;
}

void desiprostats(int num_atoms, IPRO_STATS *restrict ipro_st) {
    for (int i = 0; i < num_atoms; i++) {
        free(ipro_st->dists[i]);
    }
    free(ipro_st->dists);
}
