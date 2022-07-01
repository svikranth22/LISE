#include "lise.h"

int makeiproatom(FILE *restrict file, IPRO_ATOM *restrict atom) {
    char str1[4];
    char str2[4];
    if ((fscanf(file, "%s %d %s %f %f %f %f %d",
        str1, &atom->rno, str2, &atom->x, &atom->y, &atom->z, &atom->r, &atom->at)) == 8) {
        
        atom->aa = malloc(strlen(str1)+1);
        strcpy(atom->aa, str1);
        atom->aa[strlen(str1)] = '\0';
        atom->nom = malloc(strlen(str2)+1);
        strcpy(atom->nom, str2);
        atom->nom[strlen(str1)] = '\0';
        free(atom->aa);
        return 1;
    }
    return 0;
}

void desiproatom(IPRO_ATOM *restrict atom) {
    free(atom->aa);
    free(atom->nom);
}