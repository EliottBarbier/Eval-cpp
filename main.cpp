#include <iostream>

#include <fstream>

#include "Euler.h"
#include "Système.h"

int main(int arc,char** argv){

std::string argprin(argv[1]);
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

Cmat a;
std::vector<std::vector<double>> voulu;
std::vector<double> ligne1 {1,2};
std::vector<double> ligne2 {4,5};
std::vector<double> ligne3 {7,8};
voulu.push_back(ligne1);
voulu.push_back(ligne2);
voulu.push_back(ligne3);
a.init(voulu); //Taille (3,2)

a.affichage_mat("a");
a.transpose().affichage_mat("transposée de a");

Cmat b;
ligne1={2,4,6,8};
ligne2={1,3,5,7};//C'est 2 fois l'identité, de taille 2.
voulu={ligne1,ligne2};
b.init(voulu);
b.affichage_mat("b");

((a+a)*b).affichage_mat("(a+a)*b");
(a+ a*1.4).affichage_mat("a+1.4*a");

Cmat id;
id.identity(1,3);
id.affichage_mat("id (3,3)");
(a*a.transpose()-id).affichage_mat("a*aT - id");

Cmat pleine;
Cmat J;
Cmat diagsup;
Cmat nulle;
pleine.pleine_time_k(2.,{3,4});
J.Jp(2,3,{5,7});
diagsup.diag_sup(3,{5,4});
nulle.mat_nulle({7,9});

pleine.affichage_mat("Pleine taille 3,4 valeur 2");
J.affichage_mat("Jp taille 5,7 valeur 3");
diagsup.affichage_mat("diag sup avec valeur 3 taille 5,4");
nulle.affichage_mat("Nulle taille 7,9");

a.change_value(2,1,99);
a.affichage_mat("Valeur en bas à gauche changée en 99");

Cmat colonne;
std::vector<std::vector<double>> voulu_col;
for(int k=0;k<=8;k++){
    std::vector<double> ligne_col;
    ligne_col.push_back(1);
    voulu_col.push_back(ligne_col);
}
colonne.init(voulu_col);
colonne.affichage_mat("Colonne de 1, taille 9,1");
colonne.transpose().affichage_mat("La transposée");
(colonne.transpose()*colonne).affichage_mat("col*colT");
std::cout<< "Norme de la colonne : " << colonne.norme() <<std::endl;



    return(0);
}
else if(argprin=="test_euler"){
    


}
else if(argprin=="test_syst"){

Cmat M2;
Cmat id2;
Cmat mat_sup2;

id2.identity(-2,3); //Le -2 vient du fait que l'on a pris un certain D(x)=1
mat_sup2.diag_sup(1,{3,3});
M2 = (id2 + mat_sup2 + mat_sup2.transpose());
M2.affichage_mat("M2");

Cmat b2;
std::vector<std::vector<double>> voulub;
for(int i = 0; i<=2;i++){
    std::vector<double> ligne;
    ligne.push_back(i+1);
    voulub.push_back(ligne);
}
b2.init(voulub);
b2.affichage_mat("Ce qu'on a à droite dans le système à 3 inconnues");

Cmat nullity;
nullity.mat_nulle({3,1});
((grad_conju(M2*(-1),nullity,b2))*(-1)).affichage_mat("TEST GRAD CONJU, on devrait trouver -2.5, -4, -3.5"); //Résoud Ax=b en partant de x0, où A symétrique def positive.

}
else if(argprin=="test_lim"){
    std::cout << "deltat : " << deltat << "et deltax : "<<deltax<<std::endl;
    std::cout << "deltax^2/2 : " << pow(deltax,2)/2 <<std::endl;
    if(deltat <=  pow(deltax,2)/2){
        std::cout <<"OK" <<std::endl;
    }
    else{
        std::cout<<"Problème"<<std::endl;
    }

}
else{





//Définissons K :
Cmat K;
Cmat id;
Cmat mat_sup;
id.identity(-2,Nx); //Le -2 vient du fait que l'on a pris un certain D(x)=1
mat_sup.diag_sup(1,{Nx,Nx});
K = (id + mat_sup + mat_sup.transpose())*(1./pow(deltax,2));
//

//Définissons T0:
Cmat T0;
std::vector<std::vector<double>> voulu2;
for(int i = 0; i<=Nx-1;i++){
    std::vector<double> ligne;
    ligne.push_back(1/2 + std::sin(2*M_PI*(x_min+i*deltax)) -1/2*std::cos(2*M_PI*(x_min+i*deltax)) );
    voulu2.push_back(ligne);
}
T0.init(voulu2);
std::cout<< T0.get_shape().first << " Et " << T0.get_shape().second <<std::endl;
//

K.affichage_mat("Test K");

std::pair<std::vector<Cmat>, std::vector<double>> Euler;
Euler = Euler_explicit(T0,K,Nt,deltat,t_min,CL);
std::vector<Cmat> T;
std::vector<double> temps;
T=Euler.first;
temps=Euler.second;
T.back().affichage_mat("Retour Euler");



std::ofstream Monflux("Euler_explicite.txt");
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

std::pair<std::vector<Cmat>, std::vector<double>> Euler_imp;
Euler_imp = Euler_implicit(T0,K,Nt,deltat,t_min,CL);
std::vector<Cmat> T_imp;
std::vector<double> temps_imp;
T_imp=Euler_imp.first;
temps_imp=Euler_imp.second;
T_imp.back().affichage_mat("Retour Euler_imp");

std::ofstream Monflux_imp("Euler_implicite.txt");
if(Monflux_imp){
    int longueur = T_imp.size();
    for(int i=0; i<=T_imp.back().get_shape().first;i++){
        if (i==0){
            Monflux_imp<< "#Temps ";
        }
        else{
            Monflux_imp<<" x_"<< std::to_string(i-1);
        }
    }
    Monflux_imp<< std::endl;

    for(int k=0; k<=temps_imp.size()-1;k++ ){
        Monflux_imp<< temps_imp[k] <<" ";

        for(int i=0;i<=T_imp[k].get_shape().first-1;i++){
            Monflux_imp << T_imp[k].get_val(i,0)<< " ";
        }
        Monflux_imp<< std::endl;
        
    }
}
else{
    std::cout<<"Erreur, impossible d'ouvrir le fichier" <<std::endl;
}
}

}