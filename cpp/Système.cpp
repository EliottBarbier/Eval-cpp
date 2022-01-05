#include "Système.h"

//Cmat grad_conju(const Cmat &A, const Cmat &x0, const Cmat &b){ //A doit être r&elle symétriue positive ! Ax=b

//Cmat r = b - (A*x0);
//Cmat p=r;
//int k=0;
//Cmat x = x0;
//while(r >= 0.01){
//    float alpha = (r.transpose()*r).get_val(0,0)/((p.transpose()*A)*p).get_val(0,0);
//
//    x=x+p.scalar(alpha);
//    r=r-(A*p).scalar(alpha);
//
//    beta=(r.transpose()*r).get_val(0,0)/(r.transpose()*r).get_val(0,0);
//    p=r+p.scalar(beta);
//    k=k+1;
//}
//return(x);
//}"