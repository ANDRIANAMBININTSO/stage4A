#include <iostream>
#include <queue>
#include <utility>
#include <fstream>
#include <random>
#include <map>
#include <ostream>
#include <cassert>
#include <cmath>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Surface_mesh_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/HDVF/Hdvf_space.h>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/OSM/OSM.h>
#include <eigen3/Eigen/Dense>

//#define DEBUG

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
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits> Chain_complex;
typedef HDVF::Hdvf<Chain_complex> HDVF_type;
using PSC_flag = HDVF::PSC_flag;
using Cell_pair = HDVF::Cell_pair;

// Types Eigen
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;

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
    int distance, distance_connectedness;
    int pred;
    Operation op;
    int deg_M, deg_W, deg_MW;
    Data(int id1, int distance1, int distance_connectedness1, int pred1, Operation op1, int deg_M1, int deg_W1, int deg_MW1): id(id1), distance(distance1), distance_connectedness(distance_connectedness1), pred(pred1), op(op1), deg_M(deg_M1), deg_W(deg_W1), deg_MW(deg_MW1){}
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
    int min;
    int max;
    double mean;
    std::vector<int> hist;
    stat(int deg_min1, int deg_max1, int deg_mean1, std::vector<int> hist_M1): min(deg_min1), max(deg_max1), mean(deg_mean1), hist(hist_M1){}
    void afficher(){
        std::cout << "Min: " << min << std::endl;
        std::cout << "Max: " << max << std::endl;
        std::cout << "Mean: " << mean << std::endl;
        afficher_tab(hist);
    }
};

template <typename HdvfType>
class Hdvf_space {
public:
    typedef HdvfType Hdvf_type;
    typedef typename Hdvf_type::Coefficient_ring Coefficient_ring;
    typedef typename Hdvf_type:: Chain_complex Chain_complex;
    typedef std::map<std::vector<PSC_flag>, Data> Map_type;

protected:
    Map_type map;
    const Chain_complex& complex;
    std::string filename;
    MatrixXd shortest;
    int limit;
    Hdvf_type hdvf;
public:
    Hdvf_space(const Chain_complex& c, std::string file) : complex(c), filename(file), limit(0), hdvf(Hdvf_type(c, HDVF::OPT_FULL)) {
        // Hdvf built above
        hdvf.compute_perfect_hdvf();
        hdvf.write_hdvf_reduction("tmp/hdvf.hdvf");
//        hdvf.read_hdvf_reduction("tmp/hdvf.hdvf");

        CGAL::IO::write_VTK(hdvf, complex, "tmp/test_hdvf");

        // Run flooding
        flooding(hdvf, 1);
        std::vector<stat> vec = stat_G();
        stat M = vec[0];
        stat W = vec[1];
        stat MW = vec[2];
        stat DC = vec[3];
        stat DG = vec[4];
        stat dist = vec[5];
        stat dist_connectedness = vec[6];
        std::cout << "#######################" << filename << "#######################" << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>Stat Degree M<<<<<<<<<<<<<<<<<<" << std::endl;
        M.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat Degree W<<<<<<<<<<<<<<<<<<" << std::endl;
        W.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat Degree MW<<<<<<<<<<<<<<<<<" << std::endl;
        MW.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat Degree<<<<<<<<<<<<<<<<<" << std::endl;
        DG.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat Delta d / connectedness <<<<<<<<<<<<<<<<<" << std::endl;
        DC.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat dist (shortest path)<<<<<<<<<<<<<<<<<" << std::endl;
        dist.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat connectedness<<<<<<<<<<<<<<<<<" << std::endl;
        dist_connectedness.afficher();

        // Get HDVFs of max degree
        std::vector<size_t> max_degree_ids(get_max_degree_hdvfs(DG.max));
        std::cout << "HDVFs of max degree:" << std::endl;
        for (size_t id : max_degree_ids)
            std::cout << id << " ";
        std::cout << std::endl;
        // Get HDVFs of min degree
        std::vector<size_t> min_degree_ids(get_max_degree_hdvfs(DG.max));
        std::cout << "HDVFs of min degree:" << std::endl;
        for (size_t id : min_degree_ids)
            std::cout << id << " ";
        std::cout << std::endl;

        // Compute the shortest paths in the HDVF graph
        // Save the map
        std::map<std::vector<PSC_flag>, Data> map_flood(map);
        map.clear();
        // Init the matrix
        init_shortest();
        // Run flooding again to create the adjacency matrix
        flooding(hdvf, 1, true);
        // Compute shortest paths
        compute_shortest();
//        std::cout << "shortest:" << std::endl << shortest;
        stat_shortest();
        shortest_to_matlab();

        // Test build_hdvf_from_psc_flags
        Hdvf_type hdvf_new(build_hdvf_from_psc_flags(map.begin()->first, 1));
    }

    void init_shortest () {
        // Create a matrix of size limit x limit
        shortest = MatrixXd(limit, limit);
        for (int i=0; i < limit; ++i) {
            for (int j=0; j < limit; ++j) {
                if (i==j) shortest(i,i)=0;
                else shortest(i,j)=limit;
            }
        }
    }

    void compute_shortest () {
        bool change = true;
        int cpt=0;
        while (change) {
            change = false;
            // Pseudo product -> shortest paths (Floyd-Warshall)
            for (int k=0; k < limit; ++k) {
                for (int i=0; i < limit; ++i) {
                    for (int j=0; j < limit; ++j) {
                        int tmp(shortest(i,k)+shortest(k,j));
                        if (tmp < shortest(i,j)) {
                            shortest(i,j) = tmp;
                            change = true;
                        }
                    }
                }
            }
        }
    }

    void stat_shortest () {
        int min_shortest = limit, max_shortest = 0, somme_shortest=0;
        std::vector<int> hist_shortest(100,0);
        for (int i=0; i < limit; ++i) {
            for (int j=0; j < limit; ++j) {
                if (i != j) {
                    int coef(shortest(i,j));
                    if (coef < min_shortest)
                        min_shortest = coef;
                    if (coef > max_shortest)
                        max_shortest = coef;
                    somme_shortest += coef;
                    hist_shortest.at(coef)++;
                }
            }
        }
        stat stat_shortest(min_shortest, max_shortest, somme_shortest / (1. * limit * limit), hist_shortest);
        std::cout << "############### Stats shortests paths ###############" << std::endl;
        stat_shortest.afficher();
    }

    void shortest_to_matlab () {
        std::cout << "[";
        for (int i=0; i < limit; ++i) {
            for (int j=0; j < limit; ++j) {
                std::cout << shortest(i,j) << ", " ;
            }
            std::cout << " ; " << std::endl;
        }
        std::cout << "]";
    }

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
    std::vector<Operation> connectedness(HDVF_type& X1, HDVF_type& X_prime1, int dim){
        HDVF_type X(X1), X_prime(X_prime1);
#ifdef DEBUG
        std::cout << "---------------------DEBUT------------------" << std::endl;
#endif
        std::vector<Operation> result;
        int delta = distance(X, X_prime, dim);
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
                    Column_chain::const_iterator it = g_gamma.begin();
                    trouve = false;
                    while(it != g_gamma.end() && !(trouve)){
                        if(X_prime.psc_flag(it->first, dim) != PSC_flag::SECONDARY){
                            sigma = it->first;
                            trouve = true;
                        }
                        it++;
                    }
                    if(!trouve){
                        std::cerr << "On n'a pas trouvé de sigma" << std::endl;
                        throw(std::runtime_error("On n'a pas trouvé de sigma"));
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
                    Row_chain::const_iterator it = f_etoile_gamma.begin();
                    trouve = false;
                    while (it != f_etoile_gamma.end() && !trouve){
                        if(X_prime.psc_flag(it->first, dim) != PSC_flag::PRIMARY){
                            pi = it->first;
                            trouve = true;
                        }
                        it++;
                    }
                    if(!trouve){
                        std::cerr << "On n'a pas trouvé de pi" << std::endl;
                        throw(std::runtime_error("On n'a pas trouvé de pi"));
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
                for(int i=0; i<temp.size();i++){
                    if(X_prime.psc_flag(temp[i], dim)==PSC_flag::PRIMARY){
                        temp1.push_back(temp[i]);
                    }
                }
                int i = 0;
                trouve = false;
                while(i<temp1.size() && !trouve){
                    if(X.is_valid_pair_for_MW(pi, temp1[i], dim)){
                        sigma = temp1[i];
                        trouve = true;
                    }
                    i++;
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
                    Column_chain f_pi(CGAL::OSM::get_column(X.matrix_f(dim), pi));
                    Column_chain::const_iterator it = f_pi.begin();
                    trouve = false;
                    while(it != f_pi.end() && !trouve){
                        if(X.is_valid_pair_for_M(pi, it->first, dim)){
                            gamma = it->first;
                        }
                        it++;
                    }
                    if(!trouve){
                        std::cerr << "On n'a pas trouvé de gamma" << std::endl;
                        throw std::runtime_error("On n'a pas trouvé de gamma");
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
#ifdef DEBUG
        std::cout << "---------------------FIN------------------" << std::endl;
#endif
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


    void traiter(HDVF_type& X_origin, std::queue<HDVF_type>& a_traiter, std::map<std::vector<PSC_flag>, Data>& map,  int dim, int& id, bool compute_paths = false){
        HDVF_type X(a_traiter.front());
        bool found_M, found_W, found_MW;
        int deg_M, deg_W, deg_MW;
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
                flag_X = X.psc_flags(dim);
                flag_X1 = X1.psc_flags(dim);
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                else{
                    id++;
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                    map.insert({X1.psc_flags(dim), D});
                    a_traiter.push(X1);
#ifdef DEBUG
                    std::cout << "X " << map.at(flag_X).id << "-M(" << c.sigma << ", " << c.tau << ")-> X " << map.at(X1.psc_flags(1)).id << std::endl;
#endif
                }
                if (compute_paths) {
                    size_t idX (map.at(flag_X).id), idX1 (map.at(flag_X1).id);
                    shortest(idX, idX1) = 1;
                }
            }
        }
        if(found_W){
            for(Cell_pair c: ops_W){
                Operation o(operations::W, c.sigma, c.tau, c.dim);
                HDVF_type X1(operer(X, o));
                flag_X = X.psc_flags(dim);
                flag_X1 = X1.psc_flags(dim);
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                else{
                    id++;
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                    map.insert({X1.psc_flags(dim), D});
                    a_traiter.push(X1);
#ifdef DEBUG
                    std::cout << "X " << map.at(flag_X).id << "-W(" << c.sigma << ", " << c.tau << ")-> X " << map.at(X1.psc_flags(1)).id << std::endl;
#endif
                }
                if (compute_paths) {
                    size_t idX (map.at(flag_X).id), idX1 (map.at(flag_X1).id);
                    shortest(idX, idX1) = 1;
                }
            }
        }

        if(found_MW){
            for(Cell_pair c: ops_MW){
                Operation o(operations::MW, c.sigma, c.tau, c.dim);
                HDVF_type X1(operer(X, o));
                flag_X = X.psc_flags(dim);
                flag_X1 = X1.psc_flags(dim);
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                else{
                    id++;
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, deg_M, deg_W, deg_MW);
                    map.insert({X1.psc_flags(dim), D});
                    a_traiter.push(X1);
#ifdef DEBUG
                    std::cout << "X " << map.at(flag_X).id << "-MW(" << c.sigma << ", " << c.tau << ")-> X " << map.at(X1.psc_flags(1)).id << std::endl;
#endif
                }
                if (compute_paths) {
                    size_t idX (map.at(flag_X).id), idX1 (map.at(flag_X1).id);
                    shortest(idX, idX1) = 1;
                }
            }
        }


    }

    void flooding(HDVF_type X, int dim, bool compute_paths = false){
        std::queue<HDVF_type> a_traiter;
        int id = 0;
        Operation O0(operations::NONE, 0, 0, dim);
        bool found;
        Data D0(id, 0, 0, 0, O0, X.find_pairs_M(dim, found).size(), X.find_pairs_W(dim, found).size(), X.find_pairs_MW(dim, found).size());
        map.insert({X.psc_flags(dim),D0});
        a_traiter.push(X);
        size_t cpt(0);
        while(!(a_traiter.empty())){
#ifdef DEBUG
            std::cout << cpt++ << "===================" << std::endl;
#endif
            traiter(X, a_traiter, map, dim, id, compute_paths);
            if (!compute_paths) {
                std::vector<PSC_flag> vec = a_traiter.front().psc_flags(dim);
                map.at(vec).distance_connectedness = connectedness(X, a_traiter.front(), dim).size();
            }
            a_traiter.pop();
            if (!compute_paths)
                limit++;
        }
        std::cout << "Nombre de HDVFs: " << limit << std::endl;
    }
    void afficher_map(std::map<std::vector<PSC_flag>, Data> map){
        for(auto it=map.begin(); it!=map.end(); it++){
            it->second.afficher();
        }
    }

    std::vector<stat> stat_G(){
        std::vector<stat> result;
        int nb_elt=0, somme_dist=0, somme_dist_connectedness = 0, somme_M=0, somme_W=0, somme_MW=0, somme_dc = 0, somme_dg = 0;
        int dist_min=10000, dist_max=0, dist_mean;
        int dist_connectedness_min=10000, dist_connectedness_max=0, dist_connectedness_mean;
        int deg_min_M=10000, deg_max_M=0, deg_mean_M;
        int deg_min_W=10000, deg_max_W=0, deg_mean_W;
        int deg_min_MW=10000, deg_max_MW=0, deg_mean_MW;
        int dc_min=10000, dc_max=0, dc_mean;
        int dg_min=10000, dg_max=0, dg_mean;
        std::vector<int> hist_dist(100, 0);
        std::vector<int> hist_dist_connectedness(100, 0);
        std::vector<int> hist_M (50, 0);
        std::vector<int> hist_W (50, 0);
        std::vector<int> hist_MW (50, 0);
        std::vector<int> hist_dc(100, 0);
        std::vector<int> hist_dg(100, 0);
        for(auto it=map.begin(); it!=map.end(); it++){
            Data d = it->second;
            if (d.distance < dist_min)
                dist_min = d.distance;
            if (d.distance_connectedness < dist_connectedness_min)
                dist_connectedness_min = d.distance_connectedness;

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
            if (d.distance > dist_max)
                dist_max = d.distance;
            if (d.distance_connectedness > dist_connectedness_max)
                dist_connectedness_max = d.distance_connectedness;

            if(d.deg_M > deg_max_M){
                deg_max_M = d.deg_M;
            }
            if(d.deg_W > deg_max_W){
                deg_max_W = d.deg_W;
            }
            if(d.deg_MW > deg_max_MW){
                deg_max_MW = d.deg_MW;
            }
            /////////////////////////
            if(d.deg_MW + d.deg_M + d.deg_W > dg_max){
                dg_max = d.deg_MW + d.deg_M + d.deg_W;
            }
            if(d.deg_MW + d.deg_M + d.deg_W < dg_min){
                dg_min = d.deg_MW + d.deg_M + d.deg_W;
            }
            /////////////////////////
            if(abs(d.distance_connectedness-d.distance) < dc_min){
                dc_min = abs(d.distance_connectedness-d.distance);
            }
            if(abs(d.distance_connectedness-d.distance) > dc_max){
                dc_max = abs(d.distance_connectedness-d.distance);
            }
            hist_dist[d.distance] += 1;
            hist_dist_connectedness[d.distance_connectedness] += 1;
            hist_M[d.deg_M] += 1;
            hist_W[d.deg_W] += 1;
            hist_MW[d.deg_MW] += 1;
            hist_dg[d.deg_MW+d.deg_W+d.deg_M] += 1;
            hist_dc[abs(d.distance_connectedness-d.distance)] += 1;

            somme_dist += d.distance;
            somme_dist_connectedness += d.distance_connectedness;
            somme_M += d.deg_M;
            somme_W += d.deg_W;
            somme_MW += d.deg_MW;
            somme_dc += abs(d.distance_connectedness-d.distance);
            somme_dg += d.deg_M + d.deg_W + d.deg_MW;
            nb_elt++;
        }
        dist_mean = somme_dist / (1.0 * nb_elt);
        dist_connectedness_mean = somme_dist_connectedness / (1.0 * nb_elt);
        deg_mean_M = somme_M / (1.0 * nb_elt);
        deg_mean_W = somme_W / (1.0 * nb_elt);
        deg_mean_MW = somme_MW / (1.0 * nb_elt);
        dg_mean = somme_dg / (1.0 * nb_elt);
        dc_mean = somme_dc / (1.0 * nb_elt);
        stat stat_dist(dist_min, dist_max, dist_mean, hist_dist);
        stat stat_dist_connectedness(dist_connectedness_min, dist_connectedness_max, dist_connectedness_mean, hist_dist_connectedness);
        stat stat_M(deg_min_M, deg_max_M, deg_mean_M, hist_M);
        stat stat_W(deg_min_W, deg_max_W, deg_mean_W, hist_W);
        stat stat_MW(deg_min_MW, deg_max_MW, deg_mean_MW, hist_MW);
        stat stat_dc(dc_min, dc_max, dc_mean, hist_dc);
        stat stat_dg(dg_min, dg_max, dg_mean, hist_dg);
        result.push_back(stat_M);
        result.push_back(stat_W);
        result.push_back(stat_MW);
        result.push_back(stat_dc);
        result.push_back(stat_dg);
        result.push_back(stat_dist);
        result.push_back(stat_dist_connectedness);
        return result;
    }

    std::vector<size_t> get_max_degree_hdvfs (int max_degree) {
        std::vector<size_t> res;
        for (Map_type::const_iterator it = map.cbegin(); it != map.cend(); ++it) {
            if ((it->second.deg_M+it->second.deg_W+it->second.deg_MW) == max_degree)
                res.push_back(it->second.id);
        }
        return res;
    }

    std::vector<size_t> get_min_degree_hdvfs (int min_degree) {
        std::vector<size_t> res;
        for (Map_type::const_iterator it = map.cbegin(); it != map.cend(); ++it) {
            if ((it->second.deg_M+it->second.deg_W+it->second.deg_MW) == min_degree)
                res.push_back(it->second.id);
        }
        return res;
    }

    Hdvf_type build_hdvf_from_psc_flags (const std::vector<PSC_flag>& flags_dim, int dim) {
        std::vector<std::vector<PSC_flag> > flags(hdvf.psc_flags());
        flags[dim] = flags_dim;
        Hdvf_type hdvf_new(complex, flags, true); // complex, flags, build_reduction
        return hdvf_new;
    }
};



void compute_stat(std::string& filename, std::string nodes_file){
    HDVF::Mesh_object_io<Traits> simp;
    simp.read_simp(filename);
    simp.read_nodes_file(nodes_file);
    Chain_complex complex(simp);
    std::cout << complex << std::endl;
    Hdvf_space<HDVF_type> hs(complex, filename);
}

int main(int argc, char ** argv){

    
    std::string filename, nodes_file;
    if (argc > 3) {
        std::cerr << "usage: test_hdvf_persistence [off_file]" << std::endl;
    }
    else if (argc == 3) {
        filename = argv[1];
        nodes_file = argv[2];
        std::cout << "file: " << filename << std::endl;
        std::cout << "nodes_file: " << nodes_file << std::endl;
        compute_stat(filename, nodes_file) ;
    }
    else {
        std::vector<std::string> tab, tab_nodes;
//        tab.push_back("data/simp/three_triangles.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_d.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_1.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_inv.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_f.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_2.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

//        tab.push_back("data/simp/three_triangles_3.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");

        tab.push_back("data/simp/six_triangles.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_s.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_cs.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_cc.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_ss.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_ss_inv.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_2.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        tab.push_back("data/simp/six_triangles_f.simp");
        tab_nodes.push_back("data/simp/six_triangles.nodes");

        for (int i=0; i<tab.size(); ++i) {
            compute_stat(tab.at(i), tab_nodes.at(i));
        }
    }

    return 0;
}
