#include "test.h"

void test_matrice(){
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

a=a+a;
a.affichage_mat("a après avoir fait a=a+a");

}

void test_systeme(){
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

M2.augmente().affichage_mat("M2 augmentée");
Cmat piv = pivot_gauss(M2.augmente());
piv.affichage_mat("Matrice échelonnée réduite de la matrice augmentée (fonction pivot Gauss appliqué à la matrice augmentée)");

Cmat inv=inversion_mat(M2);
inv.affichage_mat("Inverse de M2 avec la fonction inversion_mat ");
Cmat sol_gauss = sol(M2,b2);
sol_gauss.affichage_mat("Test pivot de Gauss");

}

void test_Euler(){

}

void test_limite(const double deltat, const double deltax){
    std::cout << "deltat : " << deltat << "et deltax : "<<deltax<<std::endl;
    std::cout << "deltax^2/2 : " << pow(deltax,2)/2 <<std::endl;
    if(deltat <=  pow(deltax,2)/2){
        std::cout <<"OK" <<std::endl;
    }
    else{
        std::cout<<"Problème"<<std::endl;
    }

}