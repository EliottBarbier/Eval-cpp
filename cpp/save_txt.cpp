#include "save_txt.h"

void save_texte(const std::vector<Cmat> &T, const std::vector<double> &temps, const std::string &nom){

std::ofstream Monflux(nom);
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
}