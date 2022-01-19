#include "Cmat.h"

std::pair<std::vector<Cmat>,std::vector<double>> Euler_explicit(Cmat &T0,Cmat &K, const int &Nt, const double &deltat, const double &t0, const std::pair<double,double> &CL );
//CL sont les conditions limites aux bords à tout temps : La première en x=x_min, la deuxième en x=x_max.
std::pair<std::vector<Cmat>,std::vector<double>> Euler_implicit(Cmat &T0,Cmat &K, const int &Nt, const double &deltat, const double &t0,const std::pair<double,double> &CL );

std::pair<std::vector<Cmat>,std::vector<double>> Euler_implicit_Gauss(Cmat &T0,Cmat &K, const int &Nt, const double &deltat, const double &t0,const std::pair<double,double> &CL );

