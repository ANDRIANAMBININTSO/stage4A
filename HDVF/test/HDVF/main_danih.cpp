#include <iostream>
#include <queue>
#include <utility>
#include <fstream>
#include <random>
#include <map>
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
using HDVF_type = HDVF::Hdvf<Complex> ;
using PSC_flag = HDVF::PSC_flag;
using Cell_pair = HDVF::Cell_pair;


enum operations{
        M, W, MW, NONE
};
class Operation{
private:
    operations _op;
    size_t _sigma;
    size_t _gamma;
    int _dim;
public:
    Operation();
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
};
struct Data
{
    int id;
    int distance;
    int pred;
    Operation op;
    int deg_M, deg_W, deg_MW;
    Data(int id1, int distance1, int pred1, Operation op1, int deg_M1, int deg_W1, int deg_MW1): id(id1), distance(distance1), pred(pred1), op(op1), deg_M(deg_M1), deg_W(deg_W1), deg_MW(deg_MW1){}
    void afficher(){
        std::cout << id << ": " << distance << ", " << pred << std::endl;
    }
};
void afficher_tab(std::vector<int> tab){
    std::cout << "[ ";
    int n = tab.size();
    for(int j = 0; j<n-1; j++){
        std::cout << tab[j] << ", ";
    }
    std::cout << tab[n-1];
    std::cout << "]" << std::endl;
}
struct stat{
    int deg_min;
    int deg_max;
    double deg_mean;
    std::vector<int> hist;
    stat(int deg_min1, int deg_max1, int deg_mean1, std::vector<int> hist_M1): deg_min(deg_min1), deg_max(deg_max1), deg_mean(deg_mean1), hist(hist_M1){}
    void afficher(){
        std::cout << "Deg min: " << deg_min << std::endl;
        std::cout << "Deg max: " << deg_max << std::endl;
        std::cout << "Deg mean: " << deg_mean << std::endl;
        afficher_tab(hist);
    }
};



int distance(HDVF_type X, HDVF_type X_prime, int dim){
    int result = 0;
    std::vector<size_t> temp;
    //XC
    temp = X.psc_flags(PSC_flag::CRITICAL, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i], dim)!=PSC_flag::CRITICAL){
            result += 1;
        }
    }
    //XP
    temp = X.psc_flags(PSC_flag::PRIMARY, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i], dim)==PSC_flag::SECONDARY){
            result += 3;
        }
        else{
            if(X_prime.psc_flag(temp[i], dim)==PSC_flag::CRITICAL){
                result +=1;
            }
            
        }
    }
    //XS
    temp = X.psc_flags(PSC_flag::SECONDARY, dim);
    for (int i = 0; i<temp.size(); i++){
        if(X_prime.psc_flag(temp[i],dim)==PSC_flag::PRIMARY){
            result += 3;
        }
        else{
            if(X_prime.psc_flag(temp[i], dim)==PSC_flag::CRITICAL){
                result +=1;
            }
        }
    }

    
    return result;
}
//returns a vector of all the misaligned cells ordered by their labels.
std::vector<std::vector<size_t>> misaligned(const HDVF_type &X, const HDVF_type &X_prime, int dim){
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
void supprimer(size_t cellule, std::vector<size_t>& v){
    v.erase(std::remove(v.begin(), v.end(), cellule), v.end());
}
bool dedans(size_t sigma, Column_chain &cc){
    bool result = false;
    for(Column_chain::const_iterator it=cc.begin(); it!=cc.end(); it++){
        if(it->first==sigma){
            sigma = it->first;
            result = true;
            break;
        }
    }
    return result;
}
std::vector<Operation> connectedness(HDVF_type& X, HDVF_type& X_prime, int dim){
    std::vector<Operation> result;
    int delta = distance(X, X_prime, dim);
    std::cout << "Distance initiale: " << delta << std::endl;
    std::vector<std::vector<size_t>> misaligned_PSC = misaligned(X, X_prime, dim);
    std::vector<size_t> misaligned_P = misaligned_PSC[0], temp, temp1;
    std::vector<size_t> misaligned_S = misaligned_PSC[1];
    std::vector<size_t> misaligned_C = misaligned_PSC[2];

    bool trouve = false;
    size_t gamma, sigma, pi;
    while (delta > 0){
        if(misaligned_C.size()>0){
            gamma = misaligned_C[0];
            if(X_prime.psc_flag(gamma, dim)==PSC_flag::SECONDARY){
                const Column_chain& g_gamma(CGAL::OSM::cget_column(X.matrix_g(dim), gamma));
                for(Column_chain::const_iterator it = g_gamma.begin(); it != g_gamma.end(); ++it){
                    if(X_prime.psc_flag(it->first, dim) != PSC_flag::SECONDARY){
                        sigma = it->first;
                        break;
                    }
                }
                if(X_prime.psc_flag(sigma, dim) == PSC_flag::PRIMARY){
                    supprimer(gamma, misaligned_C);
                    supprimer(sigma, misaligned_S);
                    misaligned_C.push_back(sigma);
                    delta-=3;
                }
                else{//sigma is CRITICAL
                    supprimer(gamma, misaligned_C);
                    supprimer(sigma, misaligned_S);
                    delta-=2;
                }
                X.W(sigma, gamma, dim);
                Operation op(operations::W, sigma, gamma, dim);
                result.push_back(op);
            }
            else{
                const Row_chain& f_etoile_gamma(CGAL::OSM::cget_row(X.matrix_f(dim), gamma));
                std::cout << f_etoile_gamma << std::endl;
                for(Row_chain::const_iterator it = f_etoile_gamma.begin(); it != f_etoile_gamma.end(); ++it){
                    if(X_prime.psc_flag(it->first, dim) != PSC_flag::PRIMARY){
                        pi = it->first;
                        break;
                    }
                }
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
                X.M(pi, gamma, dim);
                Operation op(operations::M, pi, gamma, dim);
                result.push_back(op);
            }
        }
        else{
            assert(misaligned_P.size()>0);
            pi = misaligned_P[0];
            temp = X.psc_flags(PSC_flag::SECONDARY, dim);
            trouve = false;
            for(int i=0; i<temp.size();i++){
                if(X_prime.psc_flag(temp[i], dim)==PSC_flag::PRIMARY){
                    temp1.push_back(temp[i]);
                }
            }
            for(int i=0; i<temp1.size();i++){
                if(X.is_valid_pair_for_MW(pi, temp1[i], dim)){
                    sigma = temp1[i];
                    trouve = true;
                    break;
                }
            }
            if(trouve){
                X.MW(pi, sigma, dim);
                delta -= 6;
                supprimer(pi, misaligned_P);
                supprimer(sigma, misaligned_S);
                Operation op(operations::MW, pi, sigma, dim);
                result.push_back(op);
            }
            else{
                const Row_chain& f_pi(CGAL::OSM::cget_row(X.matrix_f(dim), pi));
                for(Row_chain::const_iterator it = f_pi.begin(); it != f_pi.end(); ++it){
                    if(X.is_valid_pair_for_M(it->first, gamma, dim)){
                        pi = it->first;
                        break;
                    }
                }
                X.M(pi, gamma, dim);
                delta -= 1;
                supprimer(pi, misaligned_P);
                misaligned_P.push_back(gamma);
                misaligned_C.push_back(pi);
                Operation op(operations::M, pi, gamma, dim);
                result.push_back(op);
            }
        }
    }
    return result;
}

HDVF_type operer(HDVF_type X1, Operation o){
    HDVF_type X(X1);
    if(o.operation() == operations::M){
        X.M(o.sigma(), o.gamma(), o.dimension());
    }
    else{
        if(o.operation() == operations::W){
            X.W(o.sigma(), o.gamma(), o.dimension());
        }
        else{
            if(o.operation() == operations::MW){
                X.MW(o.sigma(), o.gamma(), o.dimension());
            }
            else{}
        }
    }
    return X;
}


void traiter(std::queue<HDVF_type>& a_traiter, std::map<std::vector<PSC_flag>, Data>& map,  int dim, int& id){
    HDVF_type X(a_traiter.front());
    bool found_M, found_W, found_MW;
    int d, deg_M, deg_W, deg_MW;
    std::vector<PSC_flag> flag_X, flag_X1;
    std::vector<Cell_pair> ops_M = X.find_pairs_M(dim, found_M);
    std::vector<Cell_pair> ops_W = X.find_pairs_W(dim, found_W);
    std::vector<Cell_pair> ops_MW = X.find_pairs_MW(dim, found_MW);
    deg_M = ops_M.size();
    deg_W = ops_W.size();
    deg_MW = ops_MW.size();
    if(found_M){
        for(Cell_pair c: ops_M){
            Operation o(operations::M, c.sigma, c.tau, dim);
            HDVF_type X1(operer(X, o));
            d = distance(X1, X, dim);
            flag_X = X.psc_flags(dim);
            flag_X1 = X1.psc_flags(dim);
            if(map.find(flag_X1) != map.end()){
                if((map.at(flag_X1).distance)>(d+map.at(flag_X).distance)){
                    map.at(flag_X1).distance = d+map.at(flag_X).distance;
                    map.at(flag_X1).pred = map.at(flag_X).id;
                }
            }
            else{
                id++;
                Data D(id, d+map.at(flag_X).distance, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                map.insert({X1.psc_flags(dim), D});
                a_traiter.push(X1);
            }
        }
    }
    if(found_W){
        for(Cell_pair c: ops_W){
            Operation o(operations::W, c.sigma, c.tau, c.dim);
            HDVF_type X1(operer(X, o));
            d = distance(X1, X, dim);
            flag_X = X.psc_flags(dim);
            flag_X1 = X1.psc_flags(dim);
            if(map.find(flag_X1) != map.end()){
                if((map.at(flag_X1).distance)>(d+map.at(flag_X).distance)){
                    map.at(flag_X1).distance = d+map.at(flag_X).distance;
                    map.at(flag_X1).pred = map.at(flag_X).id;
                }
            }
            else{
                id++;
                Data D(id, d+map.at(flag_X).distance, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                map.insert({X1.psc_flags(dim), D});
                a_traiter.push(X1);
            }
        }
    }
    if(found_MW){
        for(Cell_pair c: ops_MW){
            Operation o(operations::MW, c.sigma, c.tau, c.dim);
            HDVF_type X1(operer(X, o));
            
            d = distance(X1, X, dim);
            flag_X = X.psc_flags(dim);
            flag_X1 = X1.psc_flags(dim);
            if(map.find(flag_X1) != map.end()){
                if((map.at(flag_X1).distance)>(d+map.at(flag_X).distance)){
                    map.at(flag_X1).distance = d+map.at(flag_X).distance;
                    map.at(flag_X1).pred = map.at(flag_X).id;
                }
            }
            else{
                id++;
                Data D(id, d+map.at(flag_X).distance, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                map.insert({X1.psc_flags(dim), D});
                a_traiter.push(X1);
            }
        }
    }
    a_traiter.pop();
}

std::map<std::vector<PSC_flag>, Data> flooding(HDVF_type X, int dim){
    std::queue<HDVF_type> a_traiter;
    std::map<std::vector<PSC_flag>, Data> map;
    int id = 0;
    Operation O0(operations::NONE, 0, 0, dim);
    bool found;
    Data D0(id, 0, 0, O0, X.find_pairs_M(dim, found).size(), X.find_pairs_W(dim, found).size(), X.find_pairs_MW(dim, found).size());
    map.insert({X.psc_flags(dim),D0});
    a_traiter.push(X);
    int limit = 0;
    while(!(a_traiter.empty())){
        traiter(a_traiter, map, dim, id);
        limit++;
    }
    std::cout << limit << std::endl;
    return map;
}
void afficher_map(std::map<std::vector<PSC_flag>, Data> map){
    for(auto it=map.begin(); it!=map.end(); it++){
        it->second.afficher();
    }
}

std::vector<stat> stat_G(std::map<std::vector<PSC_flag>, Data>& map){
    std::vector<stat> result;
    int nb_elt=0, somme_M=0, somme_W=0, somme_MW=0;
    int deg_min_M=10000, deg_max_M=0, deg_mean_M;
    std::map<int, int> map_M;
    int deg_min_W=10000, deg_max_W=0, deg_mean_W;
    std::map<int, int> map_W;
    int deg_min_MW=10000, deg_max_MW=0, deg_mean_MW;
    std::map<int, int> map_MW;
    std::vector<int> hist_M (500, 0);
    std::vector<int> hist_W (500, 0);
    std::vector<int> hist_MW (500, 0);
    for(auto it=map.begin(); it!=map.end(); it++){
        Data d = it->second;
        if(d.deg_M < deg_min_M){
            deg_min_M = d.deg_M;
        }
        if(d.deg_W < deg_min_W){
            deg_min_W = d.deg_W;
        }
        if(d.deg_MW < deg_min_MW){
            deg_min_MW = d.deg_MW;
        }
        ////////////////////////
        if(d.deg_M > deg_max_M){
            deg_max_M = d.deg_M;
        }
        if(d.deg_W > deg_max_W){
            deg_max_W = d.deg_W;
        }
        if(d.deg_MW > deg_max_MW){
            deg_max_MW = d.deg_MW;
        }
        hist_M[d.deg_M] += 1;
        hist_W[d.deg_W] += 1;
        hist_MW[d.deg_MW] += 1;

        somme_M += d.deg_M;
        somme_W += d.deg_W;
        somme_MW += d.deg_MW;
        nb_elt++;
    }
    deg_mean_M = somme_M / nb_elt;
    deg_mean_W = somme_W / nb_elt;
    deg_mean_MW = somme_MW / nb_elt;
    stat stat_M(deg_min_M, deg_max_M, deg_mean_M, hist_M);
    stat stat_W(deg_min_W, deg_max_W, deg_mean_W, hist_W);
    stat stat_MW(deg_min_MW, deg_max_MW, deg_mean_MW, hist_MW);
    result.push_back(stat_M);
    result.push_back(stat_W);
    result.push_back(stat_MW);
    return result;
}


int main(int argc, char ** argv){

    
    std::string chemin = "data/three_triangles.off";
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_off(chemin);

    // Build simplicial chain complex
    Complex complex(mesh);

//    // Build empty HDVF
    HDVF_type hdvf(complex, HDVF::OPT_FULL, 1);
    //HDVF_type hdvf1(complex, HDVF::OPT_FULL, 1);
    hdvf.compute_perfect_hdvf();
    //std::vector<Operation> ops = connectedness(hdvf, hdvf1, 1);

    std::map<std::vector<PSC_flag>, Data> result = flooding(hdvf, 1);
    std::vector<stat> vec = stat_G(result);
    stat M = vec[0];
    M.afficher();
    return 0;
}