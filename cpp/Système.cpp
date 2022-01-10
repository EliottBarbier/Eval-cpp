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