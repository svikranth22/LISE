#include "lise.h"


int setipro(FILE *restrict file, IPRO *restrict ipro) {
    int ret = 0;
    for (int i = 0; i < ipro->numAtoms; i++) {
        if (makeiproatom(file, ipro->atoms+ret) != 1)
            return ret;
        ret++;
    }

    return ret;
}

void createipro(int numLines, IPRO *restrict ipro) {
    ipro->numAtoms = numLines;
    ipro->atoms = (IPRO_ATOM *)malloc(numLines * sizeof(IPRO_ATOM));
}

/*
int setipro(FILE *restrict file, IPRO *restrict ipro) {
    int ret = 0;
    IPRO_ATOM atom;
    while (makeiproatom(file, &atom) == 1) {
        ret++;
        if (ret >= ipro->numAtoms) {
            ipro->atoms = realloc(ipro->atoms, ret * sizeof(atom));
            ipro->numAtoms++;
        }
        ipro->atoms[ret-1] = atom;
    }
    return ret;
}
*/
void desipro(IPRO *restrict ipro) {
    for (int i = 0; i < ipro->numAtoms; i++) {
        desiproatom(ipro->atoms+i);
    }
    free(ipro->atoms);
}