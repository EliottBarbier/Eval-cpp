#include "Euler.h"
#include "Système.h"

std::pair<std::vector<Cmat>,std::vector<double>> Euler_explicit(Cmat &T0,Cmat &K, const int &Nt, const double &deltat,const double &t0,const std::pair<double,double> &CL){
    //Cette fonction ressort un couple avec un vecteur de matrices contenant les vecteurs de températures à tout instant ainsi qu'un
    //vecteur de double contenant tout les instants correspondants.
    std::vector<double> temps;
    temps.push_back(t0);

    std::vector<Cmat> Temperature;
    Temperature.push_back(T0);
    
    Cmat T;
    T=T0;
    Cmat id;
    
    for(int i=1;i<=Nt-1;i++){
        
        id.identity(1,K.get_shape().first);
        //On applique la relation de l'énonce pour tourver T^{k+1}
        T = ((id-K*deltat)*T);

        //Les conditions limites :
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

        //On résoud le système linéaire, -M étant symétrique définie positive.
        T=grad_conju(M*(-1),T*(-1),T)*(-1); //M est symétrique définie négative c'est pour ça qu'on *(-1)
        //Les CL :
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

        
        T=sol(M,T); //SEULE DIFFERENCE AVEC l'AUTRE EULER IMPLICITE. On utilise la méthode du pivot de Gauss pour résoudre. En effet,
                    //La matrice M est toujours inversible car carré de rang n (les colonnes sont libres car échelonnées à peut prêt, en
                    // sachant que la probabilité d'avoir deux colonnes égales vaut 0 puisque on a une distribution aléatoire uniforme de flottant).
        T.change_value(0,0,CL.first);
        T.change_value(T.get_shape().first-1,0,CL.second);

        temps.push_back(i*deltat + t0);
        Temperature.push_back(T);
    }
    std::pair<std::vector<Cmat>,std::vector<double>> final{Temperature,temps};
    return(final);

}