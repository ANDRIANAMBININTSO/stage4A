#include <iostream>
#include <utility>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/HDVF/Hdvf_space.h>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/OSM/OSM.h>


namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::Zp<5, char, true> Coefficient_ring;
typedef CGAL::OSM::Sparse_chain<Coefficient_ring, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::COLUMN> Column_matrix;
typedef CGAL::OSM::Sparse_chain<Coefficient_ring, CGAL::OSM::ROW> Row_chain;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::ROW> Row_matrix;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;
using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
using HDVF_type = HDVF::Hdvf_space<Complex> ;
using PSC_flag = HDVF::PSC_flag;


enum operations{
        M, W, MW
};
class Operation{
private:
    operations _op;
    size_t _sigma;
    size_t _gamma;
    int _dim;
public:
    Operation(operations op1, size_t sigma1, size_t tau1, int dim1){
        _op = op1;
        _sigma = sigma1;
        _gamma = tau1;
        _dim = dim1;
    }
    operations operation(){return _op;}
    size_t sigma(){return _sigma;}
    size_t gamma(){return _gamma;}
    int dimension(){return _dim;}
    void display(){
        std::cout << _op << ": " << _sigma << ", " << _gamma << std::endl;
    }
};



void afficher(std::vector<size_t> tab){
    std::cout << "[ ";
    for(int j = 0; j<tab.size(); j++){
        std::cout << tab[j] << " ";
    }
    std::cout << "]" << std::endl;
}


int distance(HDVF_type X, HDVF_type X_prime, int dim){
    int result = 0;
    std::vector<size_t> temp;
    //XC
    temp = X.psc_flags(PSC_flag::CRITICAL, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(dim, i)!=PSC_flag::CRITICAL){
            result += 1;
        }
    }
    //XP
    temp = X.psc_flags(PSC_flag::PRIMARY, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(dim, i)==PSC_flag::SECONDARY){
            result += 3;
        }
        else{
            result +=1;
        }
    }
    //XS
    temp = X.psc_flags(PSC_flag::CRITICAL, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(dim, i)==PSC_flag::PRIMARY){
            result += 3;
        }
        else{
            result += 1;
        }
    }

    
    return result;
}


//returns a pair where: first is the dimension, second is the index.
std::vector<std::vector<size_t>> misaligned(HDVF_type X, HDVF_type X_prime, int dim){
    std::vector<std::vector<size_t>> result;
    std::vector<size_t> misaligned_P;
    std::vector<size_t> misaligned_S;
    std::vector<size_t> misaligned_C;
    std::vector<size_t> temp;

    //XP
    temp = X.psc_flags(PSC_flag::PRIMARY, dim);
    for(int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i], dim)!=PSC_flag::PRIMARY){
            misaligned_P.push_back(temp[i]);
        }
    }
    //XS
    temp = X.psc_flags(PSC_flag::SECONDARY, dim);
    for(int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i], dim)!=PSC_flag::SECONDARY){
            misaligned_S.push_back(temp[i]);
        }
    }
    //XC
    temp = X.psc_flags(PSC_flag::CRITICAL, dim);
    for(int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i], dim)!=PSC_flag::CRITICAL){
            misaligned_C.push_back(temp[i]);
        }
    }

    result.push_back(misaligned_P);
    result.push_back(misaligned_S);
    result.push_back(misaligned_C);

    return result;

}

void supprimer(int cellule, std::vector<size_t>& misaligned_flag){
    for(int i=0; i<(misaligned_flag.size()); i++){
        if (misaligned_flag[i] == cellule){
            misaligned_flag.erase(misaligned_flag.begin()+i);
        }
    }
}

bool dedans(size_t sigma, Column_chain cc){
    bool result = false;
    for(Column_chain::const_iterator it=cc.begin(); it!=cc.end(); it++){
        if(it->first==sigma){
            sigma = it->first;
            result = true;
        }
    }
    return result;
}

std::vector<Operation> connectedness(HDVF_type X, HDVF_type X_prime, int dim){
    std::vector<Operation> result;
    int delta = distance(X, X_prime, dim);

    std::vector<std::vector<size_t>> misaligned_PSC = misaligned(X, X_prime, dim);
    std::vector<size_t> misaligned_P = misaligned_PSC[0];
    std::vector<size_t> misaligned_S = misaligned_PSC[1];
    std::vector<size_t> misaligned_C = misaligned_PSC[2];

    int gamma, sigma, pi;
    while (delta > 0){
        if(misaligned_C.size()!=0){
            gamma = misaligned_C[0];
            if(X_prime.psc_flag(gamma, dim)==PSC_flag::SECONDARY){
                const Column_chain& g_gamma(CGAL::OSM::cget_column(X.matrix_g(dim), gamma));
                for(Column_chain::const_iterator it = g_gamma.begin(); it != g_gamma.end(); ++it){
                    sigma = it->first;
                    if(X_prime.psc_flag(sigma, dim) != PSC_flag::SECONDARY){
                        break;
                    }
                }
                X.W(sigma, gamma, dim);
                if(X_prime.psc_flag(sigma, dim) == PSC_flag::PRIMARY){
                    supprimer(gamma, misaligned_C);
                    misaligned_C.push_back(sigma);
                    supprimer(sigma, misaligned_S);
                    delta-=3;
                }
                else{
                    supprimer(gamma, misaligned_C);
                    supprimer(sigma, misaligned_S);
                    delta-=2;
                }
                    
            }
            else{
                const Row_chain& f_etoile_gamma(CGAL::OSM::cget_row(X.matrix_f(dim), gamma));
                for(Row_chain::const_iterator it = f_etoile_gamma.begin(); it != f_etoile_gamma.end(); ++it){
                    pi = it->first;
                    if(X_prime.psc_flag(pi, dim) != PSC_flag::PRIMARY){
                        break;
                    }
                }
                X.M(pi, gamma, dim);
                if(X_prime.psc_flag(pi, dim) == PSC_flag::SECONDARY){
                    supprimer(gamma, misaligned_C);
                    supprimer(pi, misaligned_P);
                    misaligned_C.push_back(pi);
                    delta -= 3;
                }
                else{
                    supprimer(gamma, misaligned_C);
                    supprimer(pi, misaligned_P);
                    delta-=2;
                }
            }
        }
        else{
            pi = misaligned_P[0];
            Column_chain dh_pi = X.htdt(pi, dim);
            for(Column_chain::const_iterator it=dh_pi.begin(); it!=dh_pi.end(); it++){
                if(X.psc_flag(it->first, dim)==PSC_flag::SECONDARY && X_prime.psc_flag(it->first, dim)==PSC_flag::PRIMARY){
                    sigma = it->first;
                    break;
                }
            }
            Column_chain hd_pi = X.hd(pi, dim);
            if(dedans(sigma, hd_pi)){
                X.MW(pi, sigma, dim);
                supprimer(pi, misaligned_P);
                supprimer(sigma, misaligned_S);
                delta -= 6;
            }
            else{
                const Row_chain& f_pi(CGAL::OSM::cget_row(X.matrix_f(dim), pi));
                gamma = f_pi.begin()->first;
                X.M(pi, gamma, dim);
                delta -= 1;
            }
        }
    }
    return result;
}

int main(int argc, char ** argv){
    
    std::string chemin = "data/three_triangles.off";
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_off(chemin);

    // Build simplicial chain complex
    Complex complex(mesh);

//    // Build empty HDVF
    HDVF_type hdvf(complex, HDVF::OPT_FULL, 1) ;
    HDVF_type hdvf1(complex, HDVF::OPT_FULL, 1) ;
    // Compute a perfect HDVF
    hdvf.compute_perfect_hdvf();
    hdvf1.compute_rand_perfect_hdvf();

    
    return 0;
}