#include <vector>
#include <utility>
#include <iostream>
#include <string>


class Cmat
{
private:
    std::pair<int,int> _taille={0,0};
    std::vector<std::vector<float>> _matrice;

public:
    Cmat(/* args */);
    ~Cmat();
public:

    void identity(const float &k,const int &taille);
    void pleine_time_k(const float &k,std::pair<int,int> &taille);
    void init(const std::vector<std::vector<float>> &voulu);
    void Jp(const float &k,const int &p, std::pair<int,int> &taille);
    void mat_nulle(const std::pair<int,int> &taille);
    void diag_sup(const float &k,const std::pair<int,int> &taille);


    std::pair<int,int> get_shape();
    
    Cmat operator+(Cmat &mat);
    Cmat operator*(Cmat &mat);
    Cmat operator-(Cmat &mat);
    //Cmat operator=(Cmat &mat); //En fait elle marche déjà...

    Cmat scalar(const float &k);
    Cmat transpose();

    void affichage_mat(const std::string &Indication);
    float get_val(const int &i,const int &j);


};

