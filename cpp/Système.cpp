#include "Système.h"
#include "Cmat.h"

Cmat grad_conju(const Cmat &A,const Cmat &x0,const Cmat &b){ //A doit être r&elle symétriue positive ! Ax=b

Cmat r = b - (A*x0); //J'ai l'impression qu'il respecte les priorités ??
Cmat p=r;

double beta;
double alpha;
int k=0;
Cmat x = x0;
Cmat r2;
while(r.norme() >= 0.0001 or k==0){
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

Cmat pivot_gauss(Cmat A){ //Pour trouver sa forme échelonnée réduite
    int r=0;
    std::pair<int,int> taille=A.get_shape();
    for(int j=0;j<=taille.second-1;j++){
        std::pair<int,int> result=max_col(r+1,j,A);
        int k = result.first;
        if(result.second != 0){
            r=r+1;
            A.change_line(k,A.get_line(k)*(1/result.second));
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
        }
    }
    return(A);
}

std::pair<int,float> max_col(const int &i_deb,const int &j, const Cmat &A){
    std::pair<int,int> taille=A.get_shape();
    double max=abs(A.get_val(0,0));
    int k=0;
    for(int i=i_deb; i<=taille.first-1;i++){
        if(abs(A.get_val(i,j)) >= max){
            max=abs(A.get_val(i,j));
            k=i;
        }
    }
    return {k,max};
}


Cmat inversion_mat(const Cmat &A){
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
    Cmat inv = inversion_mat(A);
    return(inv*b);
}