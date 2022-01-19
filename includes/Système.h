#include "Cmat.h"

Cmat grad_conju(const Cmat &A,const Cmat &x0,const Cmat &b);

Cmat pivot_gauss(Cmat A); //On change la copie de A, pas directement la référence au cas ou on en ait besoin.

std::pair<int,float> max_col(const int &i_deb, const int &j, const Cmat &A);

Cmat inversion_mat(const Cmat &A); //Inverse la matrice inversible grâce au pivot de Gauss (matrice augmentée)

Cmat sol(const Cmat &A, const Cmat&b); //La matrice K est inversible car de rang n !, donc on peut l'inverser