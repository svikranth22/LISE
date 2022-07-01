#include "lise.h"

int makeiligatom(FILE *restrict file, ILIG_ATOM *restrict atom) {
    if ((fscanf(file, "%*s %*d %*s %f %f %f",
        &atom->x, &atom->y, &atom->z)) == 3) {
        return 1;
    }
    return 0;
}