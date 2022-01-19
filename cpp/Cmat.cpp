#include "Cmat.h"

Cmat::Cmat(/* args */)
{
}
Cmat::~Cmat()
{
}

void Cmat::init(const std::vector<std::vector<double>> &voulu){
    _matrice=voulu;
    _taille={voulu.size(),voulu[0].size()};

}

void Cmat::mat_nulle(const std::pair<int,int> &taille){ //On part d'une Cmat juste initialisée.
    for(int i=0;i<=taille.first-1;i++){
        std::vector<double> ligne;
        for(int j=0;j<=taille.second-1;j++){
            ligne.push_back(0);
        }
        _matrice.push_back(ligne);
    }
    _taille=taille;
}


void Cmat::identity(const double &k,const int &taille){
    _taille={taille,taille};
    for(int i=0; i<=taille-1;i++){
        std::vector<double> lig;

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

void Cmat::pleine_time_k(const double &k,const std::pair<int,int> &taille){
    _taille={taille.first,taille.second};
    for(int i=0; i<=taille.first-1;i++){
        std::vector<double> lig(taille.second,1*k);
        _matrice.push_back(lig);
    }
}

void Cmat::Jp(const double &k, const int &p,const std::pair<int,int> &taille){ //On part du principe que p est conforme
    _taille={taille.first,taille.second};
    for(int i=0; i<=taille.first-1;i++){
        std::vector<double> lig;

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

void Cmat::diag_sup(const double &k,const std::pair<int,int> &taille){ //On part matrice juste construite vide.
    
    for(int i=0;i<=taille.first-1;i++){
        std::vector<double> ligne;
        for(int j=0;j<=taille.second-1;j++){
            if(j==i+1){
                ligne.push_back(1*k);
            }
            else{
                ligne.push_back(0);
            }
        }
        _matrice.push_back(ligne);
    }
    _taille=taille;
}


std::pair<int,int> Cmat::get_shape() const{

    return(_taille);
}

Cmat Cmat::operator+(const Cmat &mat) const{ //On part du principe qu'elles sont de même taille, pas en place.
    if(_taille.first != mat._taille.first or _taille.second != mat._taille.second ){
        throw std::length_error("Les deux matrices n'ont pas la même taille pour aditionner");
    }
    Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<double> lig;
        
        for(int j=0; j<=_taille.second-1;j++){
            lig.push_back(_matrice[i][j] + mat._matrice[i][j]);
        }
        final._matrice.push_back(lig);

    }
    final._taille=_taille;
    return(final);
}

Cmat Cmat::operator-(const Cmat &mat) const{
    if(_taille.first != mat._taille.first or _taille.second != mat._taille.second ){
        throw std::length_error("Les deux matrices n'ont pas la même taille pour soustraire");
    }
        Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<double> lig;
        
        for(int j=0; j<=_taille.second-1;j++){
            lig.push_back(_matrice[i][j] - mat._matrice[i][j]);
        }
        final._matrice.push_back(lig);
    }
    final._taille=_taille;
    return(final);
}

Cmat Cmat::operator*(const Cmat &mat) const{ //C'est bien dans l'ordre, on suppose que le nb de colonnes première correspond nb lignes deuxième.
        if(_taille.second != mat._taille.first){
        throw std::length_error("Les tailles ne correspondent pas pour la multiplication matricielle");
    }
    Cmat final;
    for(int i=0; i <= _taille.first-1;i++){

        std::vector<double> lig;
        
        for(int j=0; j<=mat._taille.second-1;j++){
            double S=0;
            for(int k=0;k<=_taille.second-1;k++){
                S=S+_matrice[i][k]*mat._matrice[k][j];
            }
            lig.push_back(S);
        }

        final._matrice.push_back(lig);

    }
    final._taille={_taille.first,mat._taille.second};
    return(final);
}

 Cmat Cmat::operator*(const double &scalar) const{
     Cmat final;
     for(std::vector<double> ligne : _matrice ){
         std::vector<double> ligne_f;
         for(double val:ligne){
            ligne_f.push_back(val*scalar);
         }
         final._matrice.push_back(ligne_f);
     }
     final._taille=_taille;
     return(final);
 }
//Cmat Cmat::operator=(Cmat &mat){
//    Cmat final;
//    final._matrice=mat._matrice;
//    final._taille=mat._taille;
//    return(final);
//}

Cmat Cmat::scalar(const double &k){
    Cmat final;
    for(int i=0; i<=_taille.first-1;i++){
        std::vector<double> ligne;
        for(int j=0; j<=_taille.second-1;j++){
            ligne.push_back(k*_matrice[i][j]);
        }
    final._matrice.push_back(ligne);
    }

    final._taille=_taille;
    return(final);
}


Cmat Cmat::transpose(){
    Cmat final;
    final.mat_nulle({_taille.second,_taille.first});

    for(int i=0;i<= final._taille.first-1;i++){
        for(int j=0;j<=final._taille.second-1;j++){
            final._matrice[i][j]=_matrice[j][i];
        }
    }
    return(final);
}



double Cmat::get_val(const int &i,const int &j) const{
    if(i>=_taille.first or j>=_taille.second){
        throw std::out_of_range("Les indices sont en dehors de la taille de la matrice");
    }
    return(_matrice[i][j]);
}



void Cmat::affichage_mat(const std::string &Indication){   
    std::cout<< Indication << std::endl;
    for(std::vector<double> ligne : _matrice){
        for(double terme : ligne){
            std::cout << terme << "|";
        }
        std::cout << std::endl;
    }
}

double Cmat::norme(){
    if(_taille.second>=2){
        throw std::length_error("La matrice n'est pas un vecteur colonne");
    }
    double S=0;
    for(int k=0;k<=_taille.first-1;k++){
        S=S+_matrice[k][0]*_matrice[k][0];
    }
    return(sqrt(S));
}

void Cmat::change_value(const int &i,const int &j,const double &value){
    if(i>=_taille.first or j>=_taille.second){
        throw std::out_of_range("Les indices sont en dehors de la taille de la matrice");
    }
    _matrice[i][j]=value;
}


Cmat Cmat::get_line(const int &i) const{
    if(i<0 or i>=_taille.first){
        throw std::out_of_range("Indice qui sort de la dimension pour get_line");
    }
    std::vector<double> ligne;
    for(int j=0; j<=_taille.second-1; j++){
        ligne.push_back( _matrice[i][j] );
    }
    Cmat final;
    (final._matrice).push_back(ligne);
    final._taille={1,_taille.second};
    return final;

}

void Cmat::change_line(const int &i, const Cmat &ligne){
    if( ligne.get_shape().first !=1 or ligne.get_shape().second != _taille.second or i<0 or i>=_taille.first){
        throw std::out_of_range("La ligne n'est pas une ligne, ou a trop de colonne, ou l'indice de la ligne a modifier sort de la matrice");
    }
    _matrice[i]=ligne._matrice[0];

}

Cmat Cmat::augmente() const{
    Cmat final;
    final._taille={_taille.first,2*_taille.second};
    
    for(int i=0; i<=_taille.first-1;i++){
        std::vector<double> ligne;
        for(int j=0;j<=2*_taille.second-1;j++){
            if(j<=_taille.second-1){
                ligne.push_back(_matrice[i][j]);
            }
            else{
                if(j==i+_taille.second){
                    ligne.push_back(1);
                }
                else{
                    ligne.push_back(0);
                }
            }
        }
        final._matrice.push_back(ligne);
    }
    return(final);
}
