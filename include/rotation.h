/*
** This is a collection of functions to manipulate matrices and triplets,
** and to calculate rotation matrices and Euler angles.
**
** Copyright (c) 2004-2010 Alexei Podtelezhnikov
*/

#ifndef ROTATION_H_
#define ROTATION_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
********** Data types
*/

typedef double matrix[3][3];
typedef vector triplet[3];

/* A trivially-copyable box holding one triplet/matrix's worth of storage, so
   it can live in a std::vector while still implicitly decaying into any
   function taking a triplet or matrix parameter -- the same underlying type
   (double[3][3]), confirmed identical, so one box covers both. No pointer
   arithmetic exists on the raw triplet* arrays this replaces (Chain::xaa/
   xaa_prev, Chaint::xaat/xaat_prev), only [i][j][k] indexing, so a box with
   just operator[] and a pointer-conversion operator is enough -- no
   constructors/destructor are declared, so it stays trivially copyable and
   Chain's implicit bitwise move stays correct. */
struct TripletBox {
	double data[3][3];
	typedef double Row[3];
	Row& operator[](int i) { return data[i]; }
	const Row& operator[](int i) const { return data[i]; }
	operator Row*() { return data; }
	operator const Row*() const { return data; }
};

/*
** Triplet and matrix manipulations
*/

void casttriplet(triplet, triplet);
void transset(matrix, triplet);
void matrixvector(vector, matrix, vector);
void vectortriplet(vector, vector, triplet);
void rotation(triplet, matrix, triplet);
void fixtriplet(triplet);
void printout(triplet);
void randvector(vector);

/*
** Rotation matrix and Euler angles
*/

double *sphereframe(vector, triplet, double, double, double);
void rotmatrix(matrix, vector, double);
void eulerset(triplet, double, double, double);
double euler_bend(triplet, triplet);
double euler_twist(triplet, triplet);
double euler_alpha(triplet, triplet);
double euler_beta(triplet, triplet);
double euler_gamma(triplet, triplet);
void tripletcmp(double *, double *, triplet, triplet);

/*
** Complex number multiplication with phase accumulation (2D-rotation)
*/

struct phasor {
	double y, x;
	int k;
};

struct phasor phasiply(struct phasor, struct phasor);
int rephase(struct phasor *);
double phase(struct phasor);


#ifdef __cplusplus
}
#endif

#endif /* ROTATION_H_ */
