#include "Cmat.h"

Cmat::Cmat(/* args */)
{
}
Cmat::~Cmat()
{
}

void Cmat::init(const std::vector<std::vector<float>> &voulu){
    _matrice=voulu;
    _taille={voulu.size(),voulu[0].size()};

}

void Cmat::identity(const float &k,const int &taille){
    _taille={taille,taille};
    for(int i=0; i<=taille-1;i++){
        std::vector<float> lig;

        for(int j=0;j<=taille-1;j++){
            if (i==j){
                lig.push_back(k*1);

            }
            else{
                lig.push_back(0);
            }
        }
        _matrice.push_back(lig);
    }
}

void Cmat::pleine_time_k(const float &k,std::pair<int,int> &taille){
    _taille={taille.first,taille.second};
    for(int i=0; i<=taille.first-1;i++){
        std::vector<float> lig(taille.second,1*k);
        _matrice.push_back(lig);
    }
}

void Cmat::Jp(const float &k, const int &p,std::pair<int,int> &taille){ //On part du principe que p est conforme
    _taille={taille.first,taille.second};
    for(int i=0; i<=taille.first-1;i++){
        std::vector<float> lig;

        for(int j=0;j<=taille.second-1;j++){
            if (i==j and j<=p-1){
                lig.push_back(k*1);

            }
            else{
                lig.push_back(0);
            }
        }
        _matrice.push_back(lig);
    }



}


std::pair<int,int> Cmat::get_shape(){

    return(_taille);
}

Cmat Cmat::operator+(Cmat &mat){ //On part du principe qu'elles sont de même taille
    Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<float> lig;
        
        for(int j=0; j<=_taille.second-1;j++){
            lig.push_back(_matrice[i][j] + mat._matrice[i][j]);
        }
        final._matrice.push_back(lig);

    }
    final._taille=_taille;
    return(final);
}

Cmat Cmat::operator-(Cmat &mat){
        Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<float> lig;
        
        for(int j=0; j<=_taille.second-1;j++){
            lig.push_back(_matrice[i][j] - mat._matrice[i][j]);
        }
        final._matrice.push_back(lig);
    }
    final._taille=_taille;
    return(final);
}

Cmat Cmat::operator*(Cmat &mat){ //C'est bien dans l'ordre, on suppose que le nb de colonnes première correspond nb lignes deuxième.
        Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<float> lig;
        
        for(int j=0; j<=_taille.second-1;j++){
            float S=0;
            for(int k=0;k<=_taille.second-1;k++){
                S=S+_matrice[i][k]*mat._matrice[k][j];
            }
            lig.push_back(S);
        }

        final._matrice.push_back(lig);

    }
    final._taille=_taille;
    return(final);
}




float Cmat::get_val(const int &i,const int &j){
    return(_matrice[i][j]);
}



void Cmat::affichage_mat(const std::string &Indication){   
    std::cout<< Indication << std::endl;
    for(std::vector<float> ligne : _matrice){
        for(float terme : ligne){
            std::cout << terme << "|";
        }
        std::cout << std::endl;
    }
}
