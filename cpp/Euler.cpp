#include "Euler.h"
#include "Système.h"

std::pair<std::vector<Cmat>,std::vector<double>> Euler_explicit(Cmat &T0,Cmat &K, const int &Nt, const double &deltat,const double &t0,const std::pair<double,double> &CL){
    std::vector<double> temps;
    temps.push_back(t0);

    std::vector<Cmat> Temperature;
    Temperature.push_back(T0);
    
    Cmat T;
    T=T0;
    Cmat id;
    
    for(int i=1;i<=Nt-1;i++){ //Il faut mettre les conditions limites ! AR car je vais pas vraiment jusqu'à x=1
        
        id.identity(1,K.get_shape().first);
        T = ((id-K*deltat)*T);

        T.change_value(0,0,CL.first);
        T.change_value(T.get_shape().first-1,0,CL.second);

        temps.push_back(i*deltat + t0);
        Temperature.push_back(T);
    }
    std::pair<std::vector<Cmat>,std::vector<double>> final{Temperature,temps};
    return(final);
}


std::pair<std::vector<Cmat>,std::vector<double>> Euler_implicit(Cmat &T0,Cmat &K, const int &Nt, const double &deltat, const double &t0,const std::pair<double,double> &CL ){
    std::vector<double> temps;
    temps.push_back(t0);

    std::vector<Cmat> Temperature;
    Temperature.push_back(T0);
    
    Cmat T;
    T=T0;

    Cmat id;
    Cmat M;
    id.identity(1.,K.get_shape().first);
    M=id+(K*deltat);
    for(int i=1;i<=Nt-1;i++){

        
        T=grad_conju(M*(-1),T*(-1),T)*(-1); //M est symétrique définie négative c'est pour ça qu'on *(-1)
        T.change_value(0,0,CL.first);
        T.change_value(T.get_shape().first-1,0,CL.second);

        temps.push_back(i*deltat + t0);
        Temperature.push_back(T);
    }
    std::pair<std::vector<Cmat>,std::vector<double>> final{Temperature,temps};
    return(final);

}

std::pair<std::vector<Cmat>,std::vector<double>> Euler_implicit_Gauss(Cmat &T0,Cmat &K, const int &Nt, const double &deltat, const double &t0,const std::pair<double,double> &CL ){
    std::vector<double> temps;
    temps.push_back(t0);

    std::vector<Cmat> Temperature;
    Temperature.push_back(T0);
    
    Cmat T;
    T=T0;

    Cmat id;
    Cmat M;
    id.identity(1.,K.get_shape().first);
    M=id+(K*deltat);
    for(int i=1;i<=Nt-1;i++){

        
        T=sol(M,T); //SEULE DIFFERENCE AVEC l'AUTRE EULER IMPLICITE.
        T.change_value(0,0,CL.first);
        T.change_value(T.get_shape().first-1,0,CL.second);

        temps.push_back(i*deltat + t0);
        Temperature.push_back(T);
    }
    std::pair<std::vector<Cmat>,std::vector<double>> final{Temperature,temps};
    return(final);

}