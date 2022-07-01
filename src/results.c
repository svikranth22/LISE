#include "lise.h"

void getresults(const IPRO_STATS *restrict ipro_st, const GRID *restrict grid, ILIG *restrict ilig, RESULT results[], int num_results) {
    float tmp;
    char tooClose = 0;
    for (int l = 0; l < num_results; l++) {
		results[l].record = 0;
		results[l].x = 0;
		results[l].y = 0;
		results[l].z = 0;

		for (int i = 0; i < ipro_st->X.range; i++) {
			for (int j = 0; j < ipro_st->Y.range; j++) {
				for (int k = 0; k < ipro_st->Z.range; k++) {
					if (grid->occ[i][j][k] == 0) {
                        if (grid->S.score[i][j][k] > results[l].record) {
                            tooClose = 0;
                            for (int m = 0; m < l; m++) {
                                tmp = sqrt(pow((float)(results[m].x - (i + ipro_st->X.min)), 2) + pow((float)(results[m].y - (j + ipro_st->Y.min)), 2) + pow((float)(results[m].z - (k + ipro_st->Z.min)), 2));
                                
                                if (tmp < bw2) {
                                    tooClose = 1;
                                    break;
                                }
                            }

                            if (tooClose == 0) {
                                results[l].record = grid->S.score[i][j][k];
                                results[l].x = i + ipro_st->X.min;
                                results[l].y = j + ipro_st->Y.min;
                                results[l].z = k + ipro_st->Z.min;
                            }
					    }
                    }
				}
			}
		}
        results[l].success = 0;
		float shortest = 100;
		if (ilig->numAtoms > 0) {
			for (int i = 0; i < ilig->numAtoms; i++) {
				tmp = sqrt(pow(ilig->atoms[i].x - (float)results[l].x, 2) + pow(ilig->atoms[i].y - (float)results[l].y, 2) + pow(ilig->atoms[i].z - (float)results[l].z, 2));
				if (shortest > tmp) {
                    shortest = tmp;
                }
			}
			if (shortest < ast) {
                results[l].success = 1;
            }
		}
        //fprintf(stdout, "HETATM    %d MN    MN     %d      %.3f  %.3f   %.3f  1.00  %.2f\n", l, l, (float)results[l].x, (float)results[l].y, (float)results[l].z, ((float)results[l].success*100));
    }
}

void printresults(const RESULT results[], int num_results, FILE *restrict results_file) {
    for (int l  = 0; l < num_results; l++) {
        fprintf(results_file, "HETATM %4d MN    MN %4d%12.3f%8.3f%8.3f  1.00  %5.2g\n", l, l, (float)results[l].x, (float)results[l].y, (float)results[l].z, ((float)results[l].success*100));
    }
}

void printdetails(const IPRO_STATS *restrict ipro_st, const GRID *restrict grid, const RESULT results[], int num_results, FILE *restrict details_file) {
    int noa = 1;
	for (int l = 0; l < 3; l++)
	{
		//cout << l << endl;
		int iu = floor(results[l].x + st) + 1, il = floor(results[l].x - st);//i:index,iu:i upper
		int ju = floor(results[l].y + st) + 1, jl = floor(results[l].y - st);//lx real coordinate lx-ipro_st->X.min translate to virtual space
		int ku = floor(results[l].z + st) + 1, kl = floor(results[l].z - st);
		if (iu > ipro_st->X.max) iu = ipro_st->X.max;
		if (ju > ipro_st->Y.max) ju = ipro_st->Y.max;
		if (ku > ipro_st->Z.max) ku = ipro_st->Z.max;
		if (il < ipro_st->X.min) il = ipro_st->X.min;
		if (jl < ipro_st->Y.min) jl = ipro_st->Y.min;
		if (kl < ipro_st->Z.min) kl = ipro_st->Z.min;
		for (int a = il; a <= iu; a++) {
			for (int b = jl; b <= ju; b++) {
				for (int c = kl; c <= ku; c++) {
					if (grid->occ[a - ipro_st->X.min][b - ipro_st->Y.min][c - ipro_st->Z.min] == 0) {
					    if (grid->loc[abs(results[l].x - a)][abs(results[l].y - b)][abs(results[l].z - c)] < st) {
						    float s = sqrt(grid->score[a - ipro_st->X.min][b - ipro_st->Y.min][c - ipro_st->Z.min] / grid->S.max * 100.0) * 10.0;
                            fprintf(details_file, "HETATM %4d  C   TYR %5d%12.3f%8.3f%8.3f  1.00 %5.2g\n", noa, l, (float)a, (float)b, (float)c, s);
						    noa++;
					    }
                    }

				}
			}
		}
	}
}