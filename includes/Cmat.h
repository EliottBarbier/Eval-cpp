#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cmath>
#include <fstream>

#pragma once //Permet de build qu'une seule fois la classe
class Cmat
{
private:
    std::pair<int,int> _taille={0,0};
    std::vector<std::vector<double>> _matrice;

public:
    Cmat(/* args */);
    ~Cmat();
public:

    //Par effets de bords :
    void identity(const double &k,const int &taille);
    void pleine_time_k(const double &k,const std::pair<int,int> &taille);
    void init(const std::vector<std::vector<double>> &voulu);
    void Jp(const double &k,const int &p, const std::pair<int,int> &taille);
    void mat_nulle(const std::pair<int,int> &taille);
    void diag_sup(const double &k,const std::pair<int,int> &taille);

    std::pair<int,int> get_shape() const;
    double get_val(const int &i,const int &j) const;
    void affichage_mat(const std::string &Indication);

    Cmat operator+(const Cmat &mat) const;
    Cmat operator*(const Cmat &mat) const;
    Cmat operator*(const double &scalar)const;
    Cmat operator-(const Cmat &mat) const;
    //Cmat operator=(Cmat &mat); //En fait elle marche déjà...

    Cmat scalar(const double &k); //Inutile maintenant
    Cmat transpose();

    Cmat get_line(const int &i) const;
    void change_line(const int &i, const Cmat &ligne); //Par effet de bord, on change A
    Cmat augmente() const;


    double norme(); //Seulement pour les vecteurs colonnes.
    void change_value(const int &i,const int &j,const double &value);

};

