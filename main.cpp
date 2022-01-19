#include <iostream>
#include <fstream>
#include <random>

#include "Euler.h"
#include "Système.h"
#include "test.h"
#include "save_txt.h"

int main(int arc,char** argv){

std::string argprin(argv[1]);
std::string argQ(argv[2]);
std::cout << argv[1] <<std::endl;

double x_max=1;
double x_min=0;
int Nx = 20; //30
int Nt=1000; //1000
double t_max=0.5;  //0.5
double t_min=0;
double deltat = (t_max-t_min)/(Nt-1);
double deltax = (x_max-x_min)/(Nx-1);
std::pair<double,double> CL{0.,0.};
//Deltat = deltax^2/2 est la limite à ne pas franchir sinon ça explose


if(argprin=="test_mat"){
    test_matrice();
    return(0);
}
else if(argprin=="test_euler"){
    test_Euler();
}
else if(argprin=="test_syst"){
    test_systeme();
}
else if(argprin=="test_lim"){
    test_limite(deltat,deltax);
}
else{


//Définissons K :
Cmat K;
if(argQ=="Q1"){
Cmat id;
Cmat mat_sup;
id.identity(-2.,Nx); //Le -2 vient du fait que l'on a pris un certain D(x)=1
mat_sup.diag_sup(1.,{Nx,Nx});
K = ((id + mat_sup + mat_sup.transpose())*(1./std::pow(deltax,2.)))*(-1); //La correction par rapport à l'énoncé consiste à
                                                                          // prendre -K. 
//
}
else if(argQ=="Q4"){
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> distr(0.5, 1.5);
    std::vector<double> D;
    for(int i=0;i<=Nx;i++){ //J'en met 1 en plus sinon on ne peut pas def K.
        D.push_back(distr(eng));
    }
    std::vector<std::vector<double>> voulu;
    for(int i=0; i<=Nx-1;i++){ 
        std::vector<double> ligne;
        for(int j=0; j<=Nx-1;j++){
            if(i==j){
                ligne.push_back(-D[i]-D[i+1]);
            }
            else if (i+1==j){
                ligne.push_back(D[i+1]);
            }
            else if(j==i-1){
                ligne.push_back(D[i]);
            }
            else{
                ligne.push_back(0.);
            }
        }
        voulu.push_back(ligne);
    }
    Cmat K2;
    K2.init(voulu);
    K = K2*(-1/pow(deltax,2));
}
//Définissons T0:
Cmat T0;
std::vector<std::vector<double>> voulu2;
for(int i = 0; i<=Nx-1;i++){
    std::vector<double> ligne;
    ligne.push_back(1/2 + std::sin(2.*M_PI*(x_min+i*deltax)) -1/2*std::cos(2.*M_PI*(x_min+i*deltax)) );
    voulu2.push_back(ligne);
}
T0.init(voulu2);
std::cout<< T0.get_shape().first << " Et " << T0.get_shape().second <<std::endl;
std::cout<< M_PI << std::endl;
//

K.affichage_mat("Test K");
std::cout<< K.get_shape().first << "Et" << K.get_shape().second <<std::endl;

//Génération des résultats pour chaque méthode : Euler explicite, Euler implicite avec et sans le gradient conjugué.

std::pair<std::vector<Cmat>, std::vector<double>> Euler= Euler_explicit(T0,K,Nt,deltat,t_min,CL);
std::vector<Cmat> T=Euler.first;
std::vector<double> temps=Euler.second;
T.back().affichage_mat("Retour Euler");
save_texte(T,temps,"Euler_explicite.txt");

if(argQ=="Q1"){
std::pair<std::vector<Cmat>, std::vector<double>> Euler_imp =Euler_implicit(T0,K,Nt,deltat,t_min,CL);
std::vector<Cmat> T_imp =Euler_imp.first;
std::vector<double> temps_imp =Euler_imp.second;
T_imp.back().affichage_mat("Retour Euler_imp");
save_texte(T_imp,temps_imp,"Euler_implicite.txt");
}

std::pair<std::vector<Cmat>, std::vector<double>> Euler_imp2 = Euler_implicit_Gauss(T0,K,Nt,deltat,t_min,CL);
std::vector<Cmat> T_imp2=Euler_imp2.first ;
std::vector<double> temps_imp2=Euler_imp2.second ;
T_imp2.back().affichage_mat("Retour Euler_imp_Gauss");
save_texte(T_imp2,temps_imp2,"Euler_implicite_Gauss.txt");




}
}