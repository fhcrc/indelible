
#include <vector>
#include <string>
#include <math.h>

extern int DiscreteGamma(std::vector<double> &, std::vector<double> &, double, double, double, int, int);
extern int DiscreteNSsites(double par[], int ngamcat, int model, std::vector<double> &output);
extern std::string newrandomtree(int ntaxa, double birth, double death, double sample, double mut, int randomseed, int option);
extern double rndgamma (double s);
extern std::vector<std::vector<double> > matexp (std::vector<std::vector<double> > Qvec, std::vector<double> &basefreqs,  double t);
extern std::vector<std::vector<double> > PMatQRev(std::vector<std::vector<double> > Qvec, std::vector<double> &basefreqs, double t);
