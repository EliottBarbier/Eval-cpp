#include "Système.h"
#include "Cmat.h"

Cmat grad_conju(const Cmat &A,const Cmat &x0,const Cmat &b){ //A doit être r&elle symétriue positive ! Ax=b

Cmat r = b - (A*x0);
Cmat p=r;

double beta;
double alpha;
int k=0;
Cmat x = x0;
Cmat r2;
while(r.norme() >= 0.0001 or k==0){ //La valeur de seuil pour r a été choisie de manière complètement arbitraire, 
                                    //pour avoir une précision suffisante sur la solution.
                   
    //L'algorithme appliqué est celui de wikipédia
    alpha = (r.transpose()*r).get_val(0,0)/((p.transpose()*A)*p).get_val(0,0);
    x=x+(p*alpha);
    r2=r-((A*p)*alpha);

    beta=(r2.transpose()*r2).get_val(0,0)/(r.transpose()*r).get_val(0,0);
    p=r2+(p*beta);
    r=r2;
    k=k+1;
}
return(x);
}

Cmat pivot_gauss(Cmat A){ //Pour trouver la forme échelonnée réduite d'une matrice A. On la recopie comme ça on peut la modifier
                          // directement sans craindre de changer la référence.
    
    //Algorithme trouvé sur Internet :
    int r=0;
    std::pair<int,int> taille=A.get_shape();
    for(int j=0;j<=taille.second-1;j++){
        if (r<=taille.first-1){ // Si r est plus petit que l'indice max de la ligne
        std::pair<int,int> result=max_col(r,j,A);
        int k = result.first;
        
        if(A.get_val(k,j) != 0){
            
            A.change_line(k,A.get_line(k)*(1/A.get_val(k,j)));
            if (k!=r){
                Cmat ligne_inter;
                ligne_inter=A.get_line(k);
                A.change_line(k,A.get_line(r));
                A.change_line(r,ligne_inter);
            }
            for(int i=0; i<=taille.first-1;i++){
                if(i!=r){
                A.change_line(i,A.get_line(i) - A.get_line(r)*(A.get_val(i,j)));
                }
            }
            r=r+1; //On incrémente r
        }
        }
    }
    return(A);
}

std::pair<int,float> max_col(const int &i_deb,const int &j, const Cmat &A){ //Permet de trouver l'indice et la valeur du maximum
//en valeur absolue d'une colonne, en partant de l'indice i_deb inclus pour les lignes.
    std::pair<int,int> taille=A.get_shape();
    double max=abs(A.get_val(i_deb,j));
    int k=i_deb;
    for(int i=i_deb; i<=taille.first-1;i++){
        if(abs(A.get_val(i,j)) >= max){
            max=abs(A.get_val(i,j));
            k=i;
        }
    }
    return {k,max};
}


Cmat inversion_mat(const Cmat &A){
    //Inverse la matrice A en utilisant la méthode de la matrice augmentée et le pivot de Gauss sur cette dernière.
    Cmat mat_augmente = A.augmente();
    Cmat piv;
    piv=pivot_gauss(mat_augmente);

    std::vector<std::vector<double>> inver_vect;
    for(int i=0; i<=piv.get_shape().first-1;i++){
        std::vector<double> ligne;
        for(int j=A.get_shape().second;j<=piv.get_shape().second -1; j++){
            ligne.push_back(piv.get_val(i,j));
        }
        inver_vect.push_back(ligne);
    }
    Cmat inv;
    inv.init(inver_vect);
    return(inv);
}

Cmat sol(const Cmat &A, const Cmat&b){
    //Résoud le système linéaire Mx=b en inversant M, puisque dans notre cas, la matrice M=Id+cst*K est toujours inversible
    Cmat inv = inversion_mat(A);
    return(inv*b);
}