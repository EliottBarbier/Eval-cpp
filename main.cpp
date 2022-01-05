#include <iostream>
#include <cmath>
#include <fstream>

#include "Euler.h"
//#include "Système.h"

int main(int arc,char** argv){

float x_max=1; //x_min=0
int Nx = 30;
int Nt=1000;
float t_max=0.5; //t_min=0
float t0=0;
float deltat = t_max/Nt;
float deltax = x_max/Nx;


//Définissons K :
Cmat K;
Cmat id;
Cmat mat_sup;

id.identity(-2,Nx); //Le -2 vient du fait que l'on a pris un certain D(x)=1
mat_sup.diag_sup(1,{Nx,Nx});
Cmat trans=mat_sup.transpose();
K = ((id + mat_sup) + trans);
//

//Définissons T0:
Cmat T0;
std::vector<std::vector<float>> voulu2;
for(int i = 0; i<=Nx-1;i++){
    std::vector<float> ligne;
    ligne.push_back(1/2 + std::sin(2*M_PI*i*deltax) -1/2*std::cos(2*M_PI*i*deltax) );
    voulu2.push_back(ligne);
}
T0.init(voulu2);
std::cout<< T0.get_shape().first << " Et " << T0.get_shape().second <<std::endl;
//

K.affichage_mat("Test K");

Cmat a2;
std::vector<std::vector<float>> voulute;
std::vector<float> ligne1te {1,2};
std::vector<float> ligne2te {4,5};
std::vector<float> ligne3te {7,8};
voulute.push_back(ligne1te);
voulute.push_back(ligne2te);
voulute.push_back(ligne3te);
a2.init(voulute);

float scalar=2;

Cmat c2;
c2=a2*scalar;
c2.affichage_mat("a*scalar");










std::pair<std::vector<Cmat>, std::vector<float>> Euler;
Euler = Euler_explicit(T0,K,Nt,deltat,t0);
std::vector<Cmat> T;
std::vector<float> temps;
T=Euler.first;
temps=Euler.second;
T.back().affichage_mat("Retour Euler");



std::ofstream Monflux("Euler.txt");
if(Monflux){
    int longueur = T.size();
    for(int i=0; i<=T.back().get_shape().first;i++){
        if (i==0){
            Monflux<< "#Temps ";
        }
        else{
            Monflux<<" x_"<< std::to_string(i-1);
        }
    }
    Monflux<< std::endl;

    for(int k=0; k<=temps.size()-1;k++ ){
        Monflux<< temps[k] <<" ";

        for(int i=0;i<=T[k].get_shape().first-1;i++){
            Monflux << T[k].get_val(i,0)<< " ";
        }
        Monflux<< std::endl;
        
    }
}
else{
    std::cout<<"Erreur, impossible d'ouvrir le fichier" <<std::endl;
}






Cmat a;
std::vector<std::vector<float>> voulu;
std::vector<float> ligne1 {1,2};
std::vector<float> ligne2 {4,5};
std::vector<float> ligne3 {7,8};
voulu.push_back(ligne1);
voulu.push_back(ligne2);
voulu.push_back(ligne3);
a.init(voulu);

Cmat b;
b.identity(2,2); //C'est 2 fois l'identité, de taille 2.

Cmat c;
c=a*b;
c.affichage_mat("a*b");


//Cmat d;
//d=a+b;
//d.affichage_mat("a+b");

Cmat e;
e=a.transpose();
e.affichage_mat("transpose de a");

Cmat f;
f.diag_sup(3, {3,3});
f.affichage_mat("diagsup");

f=e;
f.affichage_mat("test égal");
}