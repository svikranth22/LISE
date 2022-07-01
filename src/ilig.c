#include "lise.h"

int setilig(FILE *restrict file, ILIG *restrict ilig) {
    int ret = 0;
    for (int i = 0; i < ilig->numAtoms; i++) {
        if (makeiligatom(file, ilig->atoms+ret) != 1)
            return ret;
        ret++;
    }

    return ret;
}

void createilig(int numLines, ILIG *restrict ilig) {
    ilig->numAtoms = numLines;
    ilig->atoms = (ILIG_ATOM *)malloc(numLines * sizeof(ILIG_ATOM));
}
 /*
int setilig(FILE *restrict file, ILIG *restrict ilig) {
    int ret = 0;
    ILIG_ATOM atom;
    while (makeiligatom(file, &atom) == 1) {
        ret++;
        if (ret >= ilig->numAtoms) {
            ilig->atoms = realloc(ilig->atoms, ret * sizeof(atom));
            ilig->numAtoms++;
        }
        ilig->atoms[ret-1] = atom;
    }
    return ret;
}
*/
void desilig(ILIG *restrict ilig) {
    free(ilig->atoms);
}