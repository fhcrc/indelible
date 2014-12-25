#include "MersenneTwister.h"

extern int Zipf(double, double);
extern int oldZipf(double);
extern int randnegbin(int, double);
extern int oldrandnegbin(int, double);

class MTRand;   // forward declaration
extern MTRand       mtrand1;
extern MTRand       mtrand2;
extern unsigned int imax;
extern int          idum;
extern double       Zq1, Zq2, ZHx0, Zs, ZHimax;
extern double       myrand;
extern double expmyrand(double myrate);

#define newZipf(a)    Zipf(a, 1)

// Zipfian random number generation

#define H(x, q1, q2, v)     exp(q1 * log(v + x)) * q2

#define H1(x, q1, q2, v)    - v + exp(q2 * log((1 - q) * x))
