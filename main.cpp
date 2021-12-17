#include <iostream>

#include "Cmat.h"


int main(int arc,char** argv){

Cmat a;
std::vector<std::vector<float>> voulu;
std::vector<float> ligne1 {1,2};
std::vector<float> ligne2 {4,5};
voulu.push_back(ligne1);
voulu.push_back(ligne2);
a.init(voulu);

Cmat b;
b.identity(2,2); //C'est 2 fois l'identité, de taille 2.

Cmat c;
c=a*b;
c.affichage_mat("a*b");


Cmat d;
d=a+b;
d.affichage_mat("a+b");


}