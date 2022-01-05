#include "Euler.h"


std::pair<std::vector<Cmat>,std::vector<float>> Euler_explicit(Cmat &T0,Cmat &K, const int &Nt, const float &deltat,const float &t0){
    std::vector<float> temps;
    temps.push_back(t0);

    std::vector<Cmat> Temperature;
    Temperature.push_back(T0);
    
    Cmat T;
    T=T0;
    for(int i=1;i<=Nt-1;i++){
        Cmat prod=K*T;
        Cmat fin = prod.scalar(deltat);
        T = T - fin;
        temps.push_back(i*deltat + t0);
        Temperature.push_back(T);
    }
    std::pair<std::vector<Cmat>,std::vector<float>> final{Temperature,temps};
    return(final);
}
