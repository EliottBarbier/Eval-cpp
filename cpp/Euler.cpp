#include "Euler.h"


Cmat Euler_explicit(Cmat &T0,Cmat &K, const int &Nt, const float &deltat){
    Cmat T;
    T=T0;
    for(int i=1;i<=Nt-1;i++){
        Cmat prod=K*T;
        Cmat fin = prod.scalar(deltat);
        T = T - fin;
    }
    return(T);
}
