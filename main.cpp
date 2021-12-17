#include <iostream>
#include <cmath>
#include <fstream>


#include "Euler.h"

int main(int arc,char** argv){

float x_max=1;
int Nx = 100;
int Nt=1000;
float t_max=0.5;
float deltat = t_max/Nt;
float deltax = x_max/Nx;


//Définissons K :
Cmat K;
Cmat id;
Cmat mat_sup;
id.identity(-2,Nx);
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
//

//K.affichage_mat("Test K");
Cmat T;
T=Euler_explicit(T0,K,Nt,deltat);
T.affichage_mat("Retour Euler");

std::ofstream Monflux("Euler.txt");
if(Monflux){
    std::pair<int,int> taille= T.get_shape();

    for(int i=0;i<=taille.first-1;i++){
        Monflux << T.get_val(i,0) <<std::endl;
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