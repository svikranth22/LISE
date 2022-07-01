#ifndef LISE_H
#define LISE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern const float ot; //occupied grid threshold
extern const float st; //sphere threshold
extern const float ast;  //score assignment threshold
extern const float sat; //0.025;//0? //surface area threshould
extern const float Eft;
extern const float w;
extern const float bw2; //two sphere overlap threshold
extern const float MIN;
extern const float MAX;

extern float Fb[14][14][14];

extern const int rank;

//iproatom
typedef struct IPRO_ATOM {
	char* aa; //amino acid name-3 letter
	int rno; //residue serial number
	char* nom; //atom name
	float x; //x coordinates
	float y; //y coordinates
	float z; //z coordinates
	float r; //r radius
	int at; // atom type
	//float sa; //surface area
}IPRO_ATOM;

int makeiproatom(FILE *restrict, IPRO_ATOM *restrict atom);

void desiproatom(IPRO_ATOM *restrict atom);

//ipro
typedef struct IPRO {
	int numAtoms;
	IPRO_ATOM *atoms;
}IPRO;

void createipro(int numLines, IPRO *restrict ipro);

int setipro(FILE *restrict file, IPRO *restrict ipro);

void desipro(IPRO *restrict ipro);

//iligatom
typedef struct ILIG_ATOM {
	float x; //x coordinates
	float y; //y coordinates
	float z; //z coordinates
}ILIG_ATOM;

int makeiligatom(FILE *restrict, ILIG_ATOM *restrict atom);

//ilig
typedef struct ILIG {
	int numAtoms;
	ILIG_ATOM *atoms;
}ILIG;

void createilig(int numLines, ILIG *restrict ilig);

int setilig(FILE *restrict file, ILIG *restrict ilig);

void desilig(ILIG *restrict ilig);

//stats
typedef struct STATS {
	int min;
	int max;
	int range;
}STATS;

void initstats(STATS *restrict st);

void setstats(int n, STATS *restrict st);

void finstats(STATS *restrict st);

//iprostats
typedef struct IPRO_STATS {
    STATS X;
    STATS Y;
    STATS Z;
    int **dists;
}IPRO_STATS;

void createiprostats(int num_atoms, IPRO_STATS *restrict ipro_st);

int setiprostats(const IPRO *ipro, IPRO_STATS *restrict ipro_st);

void desiprostats(int num_atoms, IPRO_STATS *restrict ipro_st);

//grid
typedef struct GRID {
	char ***occ; //grid occupied
	float ***score; //grid score
	float loc[10][10][10];

	struct Sphere {
		float max;
		float ***score;
	}S;

}GRID;

void creategrid(const IPRO_STATS *restrict ipro_st, GRID *restrict grid);

int setgridocc(const IPRO *restrict ipro, const IPRO_STATS *restrict ipro_st, GRID *restrict grid);

int setgridscore(const IPRO *restrict ipro, const IPRO_STATS *restrict ipro_st, GRID *restrict grid);

void setgridloc(GRID *restrict grid);

int setgridSscore(const IPRO_STATS *restrict ipro_st, GRID *restrict grid);

void desgrid(const IPRO_STATS *restrict ipro_st, GRID *restrict grid);


typedef struct RESULT {
	float record;
	int x;
	int y;
	int z;
	char success;
} RESULT;

void getresults(const IPRO_STATS *restrict ipro_st, const GRID *restrict grid, ILIG *restrict ilig, RESULT results[], int num_results);

void printresults(const RESULT results[], int num_results, FILE *restrict results_file);

void printdetails(const IPRO_STATS *restrict ipro_st, const GRID *restrict grid, const RESULT results[], int num_results, FILE *restrict details_file);

#endif