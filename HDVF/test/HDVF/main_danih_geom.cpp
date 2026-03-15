#include <iostream>
#include <queue>
#include <utility>
#include <fstream>
#include <random>
#include <map>
#include <ostream>
#include <cassert>
#include <cmath>
#include <cstring>
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
#include "distrib.h"

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
    void afficher() const {
        if (_op == M)
            std::cout << "M(";
        else if (_op == W)
            std::cout << "W(";
        else
            std::cout << "MW(";

        std::cout << _sigma << ", " << _gamma << ")" << std::endl;
    }
};

struct Data
{
    size_t id;
    size_t distance, distance_connectedness;
    size_t pred;
    Operation op;
    size_t deg_M, deg_W, deg_MW;
    size_t len_hom, len_cohom;
    Data(size_t id1, size_t distance1, size_t distance_connectedness1, size_t pred1, Operation op1, size_t deg_M1, size_t deg_W1, size_t deg_MW1, size_t len_hom1, size_t len_cohom1): id(id1), distance(distance1), distance_connectedness(distance_connectedness1), pred(pred1), op(op1), deg_M(deg_M1), deg_W(deg_W1), deg_MW(deg_MW1), len_hom(len_hom1), len_cohom(len_cohom1) {}
    void afficher(){
        std::cout << id << ": " << distance << ", " << pred << std::endl;
    }
};

template <typename T>
void afficher_tab(std::vector<T> tab){
    std::cout << "[ ";
    int n = tab.size();
    if (n>0) {
        for(int j = 0; j<n-1; j++){
            std::cout << tab[j] << ", ";
        }
        std::cout << tab[n-1];
    }
    std::cout << "]" << std::endl;
}

template <typename T>
struct stat{
    T min;
    T max;
    double mean;
    std::vector<size_t> hist;
    std::vector<T> hist_labels;
    stat(T min1, T max1, T mean1, std::vector<size_t> hist1, std::vector<T> hist_labels1): min(min1), max(max1), mean(mean1), hist(hist1), hist_labels(hist_labels1) {}
    void afficher(){
        std::cout << "Min: " << min << std::endl;
        std::cout << "Max: " << max << std::endl;
        std::cout << "Mean: " << mean << std::endl;
        afficher_tab(hist_labels);
        afficher_tab(hist);
    }
};

template <typename HdvfType>
class Hdvf_space {
public:
    typedef HdvfType Hdvf_type;
    typedef typename Hdvf_type::Coefficient_ring Coefficient_ring;
    typedef typename Hdvf_type:: Chain_complex Chain_complex;
    typedef std::vector<PSC_flag> Flag_type;
    typedef std::map<std::vector<PSC_flag>, Data> Map_type;

protected:
    int dim;
    Map_type map;
    std::vector<Flag_type> id_to_flags;
    const Chain_complex& complex;
    std::string filename;
    MatrixXd shortest;
    int limit;
    Hdvf_type hdvf;
public:
    Hdvf_space(const Chain_complex& c, int dim1, std::string file) : complex(c), dim(dim1), filename(file), limit(0), hdvf(Hdvf_type(c, HDVF::OPT_FULL)) {
        // Hdvf built above
        hdvf.compute_perfect_hdvf();
//        hdvf.write_hdvf_reduction("tmp/hdvf.hdvf");
//        hdvf.read_hdvf_reduction("tmp/hdvf.hdvf");

        CGAL::IO::write_VTK(hdvf, complex, compute_name(file, "_init"));

        // Run flooding
        flooding(hdvf, dim);
        std::vector<stat<size_t> > vec = compute_stats();
        stat<size_t> M = vec[0];
        stat<size_t> W = vec[1];
        stat<size_t> MW = vec[2];
        stat<size_t> DC = vec[3];
        stat<size_t> DG = vec[4];
        stat<size_t> dist = vec[5];
        stat<size_t> dist_connectedness = vec[6];
        stat<size_t> len_hom = vec[7];
        stat<size_t> len_cohom = vec[8];
        stat<size_t> len_hom_cohom = vec[9];

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
        std::cout << ">>>>>>>>>>>>>>>>Stat len hom generators<<<<<<<<<<<<<<<<<" << std::endl;
        len_hom.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat len cohom generators<<<<<<<<<<<<<<<<<" << std::endl;
        len_cohom.afficher();
        std::cout << ">>>>>>>>>>>>>>>>Stat len hom+cohom generators<<<<<<<<<<<<<<<<<" << std::endl;
        len_hom_cohom.afficher();

        // Get HDVFs of max degree
        std::vector<size_t> max_degree_ids(get_max_degree_hdvfs(DG.max));
        std::cout << "HDVFs of max degree:" << std::endl;
        for (size_t id : max_degree_ids)
            std::cout << id << " ";
        std::cout << std::endl;
        // Get HDVFs of min degree
        std::vector<size_t> min_degree_ids(get_min_degree_hdvfs(DG.min));
        std::cout << "HDVFs of min degree:" << std::endl;
        for (size_t id : min_degree_ids)
            std::cout << id << " ";
        std::cout << std::endl;

        // Get HDVFs of min len generators
        std::vector<size_t> min_hom_gene_ids(get_min_hom_generators(len_hom.min));
        std::cout << "HDVFs of min hom length: " << min_hom_gene_ids.size() << " HDVFs" << std::endl;
        for (size_t id : min_hom_gene_ids)
            std::cout << id << " ";
        std::cout << std::endl;

//        for (size_t id : min_hom_gene_ids) {
//            Data tmp(map.at(id_to_flags.at(id)));
//            std::cout << "id/deg_M/deg_W/deg_MW/deg//len_hom/len_cohom: " << id << " / " << tmp.deg_M << " / " << tmp.deg_W << " / " << tmp.deg_MW << " / " << tmp.deg_M+tmp.deg_W+tmp.deg_MW << " // " << tmp.len_hom << " / " << tmp.len_cohom ;
//            std::cout << std::endl;
//        }


        // Compute the shortest paths in the HDVF graph
        // Save the map and id_to_flags
        std::map<std::vector<PSC_flag>, Data> map_flood(map);
        map.clear();
        id_to_flags.clear();
        // Init the matrix
        init_shortest();
        // Run flooding again to create the adjacency matrix
        flooding(hdvf, dim, true);
        // Compute shortest paths
        compute_shortest();
//        std::cout << "shortest:" << std::endl << shortest;
        stat_shortest();
        std::string matlab_file(compute_name(filename,".m"));
        shortest_to_matlab(matlab_file);

        // Restore the map (for connectedness computation)
        map = map_flood;

        // --------- Other computations


        // Max degrees HDVFs (N elements max)
//        const int N=1;
//        for (int i=0; (i<max_degree_ids.size()) && (i<N); ++i) {
//            Hdvf_type hdvf_new(build_hdvf_from_psc_flags(id_to_flags.at(max_degree_ids.at(i)), dim));
//            std::string suffix("_max_"+std::to_string(i)), out_hdvf_vtk_files(compute_name(filename, suffix));
//            CGAL::IO::write_VTK(hdvf_new, complex, out_hdvf_vtk_files);
//        }

//        // Min degrees HDVFs (5 elements max)
//        for (int i=0; (i<min_degree_ids.size()) && (i<5); ++i) {
//            Hdvf_type hdvf_new(build_hdvf_from_psc_flags(id_to_flags.at(min_degree_ids.at(i)), dim));
//            std::string suffix("_min_"+std::to_string(i)), out_hdvf_vtk_files(compute_name(filename, suffix));
//            CGAL::IO::write_VTK(hdvf_new, complex, out_hdvf_vtk_files);
//        }
//
//        // One min hom generator
////        Hdvf_type hdvf_new(build_hdvf_from_psc_flags(id_to_flags.at(min_hom_gene_ids.at(0)), dim));
//        Hdvf_type hdvf_new(build_hdvf_from_psc_flags(id_to_flags.at(min_hom_gene_ids.at(min_hom_gene_ids.size()-1)), dim));
//        std::string suffix("_min_hom"), out_hdvf_vtk_files(compute_name(filename, suffix));
//        CGAL::IO::write_VTK(hdvf_new, complex, out_hdvf_vtk_files);

        // Computation of a pair of paths of max delta
        std::pair<std::vector<Operation>, std::vector<Operation> > p(get_pair_paths_delta(hdvf, DC.max));
        std::cout << "shortest path length: " << p.first.size() << " - connectedness length: " << p.second.size() << std::endl;
        std::cout << "--- shortest path" << std::endl;
        export_path_to_vtk(hdvf, p.first, "tmp/hdvf_shortest_path");
        std::cout << "--- connectedness path" << std::endl;
        export_path_to_vtk(hdvf, p.second, "tmp/hdvf_connect_path");
    }

    std::string compute_name(std::string filename, std::string suffix) {
        // Cut wrt '/'
        std::stringstream str_stream1(filename);
        std::string segment;
        std::vector<std::string> seglist;

        while(std::getline(str_stream1, segment, '/')) {
           seglist.push_back(segment);
        }
        size_t n(seglist.size());
        std::string file(seglist.at(n-1)); // File name without leading path
        // Remove .xyz suffix
        std::stringstream str_stream2(file);
        while(std::getline(str_stream1, segment, '.')) {}

        return ("tmp/"+segment+suffix);
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
        distrib<size_t> dist_shortest, dist_row, dist_rows;
        for (int i=0; i < limit; ++i) {
            dist_row.clear();
            for (int j=0; j < limit; ++j) {
                if (i != j) {
                    int coef(shortest(i,j));
                    dist_shortest.add_data(coef);
                    dist_row.add_data(coef);
                }
            }
            dist_rows.add_data(dist_row.get_max());
        }

        stat<size_t> stat_shortest(dist_shortest.get_min(), dist_shortest.get_max(), dist_shortest.get_mean(), dist_shortest.get_hist(), dist_shortest.get_hist_labels());
        std::cout << "############### Stats shortests paths ###############" << std::endl;
        stat_shortest.afficher();
        std::cout << "minimum outgoing path: " << dist_rows.get_min() << std::endl;
    }

    void shortest_to_matlab (std::string filename) {
        std::ofstream out_file ( filename, std::ios::out | std::ios::trunc);
        if ( ! out_file . good () ) {
            std::cerr << "Out fatal Error:\n  " << filename << " not found.\n";
            throw std::runtime_error("File Parsing Error: File not found");
        }
        out_file << "A = [";
        for (int i=0; i < limit; ++i) {
            for (int j=0; j < limit; ++j) {
                out_file << shortest(i,j) << ", " ;
            }
            out_file << " ; " << std::endl;
        }
        out_file << "]" << std::endl;
        out_file.close();
    }

    int distance(const HDVF_type& X, const HDVF_type& X_prime, int dim){
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
    // Compute connectedness path
    std::vector<Operation> connectedness(const HDVF_type& X1, const HDVF_type& X_prime1, int dim){
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

    HDVF_type operer(const HDVF_type& X1, Operation o){
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
        // Data in the map is updated in three steps
        // 1) first meeting: give and id, set the distance, insert in the queue
        // 2) head of the queue: set degrees and visit sons
        // 3) further meetings: update the distance

        // X becomes the head of the queue -> set degrees
        HDVF_type X(a_traiter.front());
        bool found_M, found_W, found_MW;
        size_t deg_M, deg_W, deg_MW;
        std::vector<PSC_flag> flag_X, flag_X1;
        std::vector<Cell_pair> ops_M = X.find_pairs_M(dim, found_M);
        std::vector<Cell_pair> ops_W = X.find_pairs_W(dim, found_W);
        std::vector<Cell_pair> ops_MW = X.find_pairs_MW(dim, found_MW);
        deg_M = ops_M.size();
        deg_W = ops_W.size();
        deg_MW = ops_MW.size();
        flag_X = X.psc_flags(dim);
        map.at(flag_X).deg_M = deg_M;
        map.at(flag_X).deg_W = deg_W;
        map.at(flag_X).deg_MW = deg_MW;
        map.at(flag_X).len_hom = deg_W+X.number_of_cells_by_flag(HDVF::CRITICAL, dim);
        map.at(flag_X).len_cohom = deg_M+X.number_of_cells_by_flag(HDVF::CRITICAL, dim);
        // For sons
        if(found_M){
            for(Cell_pair c: ops_M){
                Operation o(operations::M, c.sigma, c.tau, dim);
                HDVF_type X1(operer(X, o));
                flag_X1 = X1.psc_flags(dim);
                // If already met: update the distance
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                // First meeting: set an id and insert in the map with degrees set to 0 (not defined)
                else{
                    id_to_flags.push_back(X1.psc_flags(dim));
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, 0, 0, 0, 0, 0);
                    id++;
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
                flag_X1 = X1.psc_flags(dim);
                // If already met: update the distance
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                // First meeting: set an id and insert in the map with degrees set to 0 (not defined)
                else{
                    id_to_flags.push_back(X1.psc_flags(dim));
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, 0, 0, 0, 0, 0);
                    id++;
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
                flag_X1 = X1.psc_flags(dim);
                // If already met: update the distance
                if(map.find(flag_X1) != map.end()){
                    if((map.at(flag_X1).distance)>(1+map.at(flag_X).distance)){
                        map.at(flag_X1).distance = 1+map.at(flag_X).distance;
                        map.at(flag_X1).pred = map.at(flag_X).id;
                    }
                }
                // First meeting: set an id and insert in the map with degrees set to 0 (not defined)
                else{
                    id_to_flags.push_back(X1.psc_flags(dim));
                    Data D(id, 1+map.at(flag_X).distance, 0, map.at(flag_X).id, o, 0, 0, 0, 0, 0);
                    id++;
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
        size_t deg_M(X.find_pairs_M(dim, found).size()), deg_W(X.find_pairs_W(dim, found).size());
        Data D0(id, 0, 0, 0, O0, deg_M, deg_W, X.find_pairs_MW(dim, found).size(), deg_W+X.number_of_cells_by_flag(HDVF::CRITICAL, dim), deg_M+X.number_of_cells_by_flag(HDVF::CRITICAL, dim));
        id_to_flags.push_back(X.psc_flags(dim));
        id++;
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

    std::vector<stat<size_t> > compute_stats(){
        distrib<size_t> distrib_dist, distrib_connectedness, distrib_M, distrib_W, distrib_MW, distrib_deg, distrib_delta, distrib_hom ,distrib_cohom, distrib_hom_cohom;

        // Visit the map and record data in various distributions
        for(auto it=map.begin(); it!=map.end(); it++){
            Data d(it->second);
            Flag_type flag(it->first);
            distrib_dist.add_data(d.distance);
            distrib_connectedness.add_data(d.distance_connectedness);
            distrib_delta.add_data(d.distance_connectedness-d.distance);
            distrib_M.add_data(d.deg_M);
            distrib_W.add_data(d.deg_W);
            distrib_MW.add_data(d.deg_MW);
            distrib_deg.add_data(d.deg_M+d.deg_W+d.deg_MW);
            distrib_hom.add_data(d.len_hom);
            distrib_cohom.add_data(d.len_cohom);
            distrib_hom_cohom.add_data(d.len_hom+d.len_cohom);
        }
// Export degree stats to matlab
        std::cout << "====> Degrees " << std::endl;
        distrib_deg.to_matlab();

        std::vector<stat<size_t> > result;
        stat<size_t> stat_dist(distrib_dist.get_min(), distrib_dist.get_max(), distrib_dist.get_mean(), distrib_dist.get_hist(), distrib_dist.get_hist_labels());
        stat<size_t> stat_connectedness(distrib_connectedness.get_min(), distrib_connectedness.get_max(), distrib_connectedness.get_mean(), distrib_connectedness.get_hist(), distrib_connectedness.get_hist_labels());
        stat<size_t> stat_M(distrib_M.get_min(), distrib_M.get_max(), distrib_M.get_mean(), distrib_M.get_hist(), distrib_M.get_hist_labels());
        stat<size_t> stat_W(distrib_W.get_min(), distrib_W.get_max(), distrib_W.get_mean(), distrib_W.get_hist(), distrib_W.get_hist_labels());
        stat<size_t> stat_MW(distrib_MW.get_min(), distrib_MW.get_max(), distrib_MW.get_mean(), distrib_MW.get_hist(), distrib_MW.get_hist_labels());
        stat<size_t> stat_deg(distrib_deg.get_min(), distrib_deg.get_max(), distrib_deg.get_mean(), distrib_deg.get_hist(), distrib_deg.get_hist_labels());
        stat<size_t> stat_delta(distrib_delta.get_min(), distrib_delta.get_max(), distrib_delta.get_mean(), distrib_delta.get_hist(), distrib_delta.get_hist_labels());
        stat<size_t> stat_hom(distrib_hom.get_min(), distrib_hom.get_max(), distrib_hom.get_mean(), distrib_hom.get_hist(), distrib_hom.get_hist_labels());
        stat<size_t>stat_cohom(distrib_cohom.get_min(), distrib_cohom.get_max(), distrib_cohom.get_mean(), distrib_cohom.get_hist(), distrib_cohom.get_hist_labels());
        stat<size_t> stat_hom_cohom(distrib_hom_cohom.get_min(), distrib_hom_cohom.get_max(), distrib_hom_cohom.get_mean(), distrib_hom_cohom.get_hist(), distrib_hom_cohom.get_hist_labels());

        result.push_back(stat_M);
        result.push_back(stat_W);
        result.push_back(stat_MW);
        result.push_back(stat_delta);
        result.push_back(stat_deg);
        result.push_back(stat_dist);
        result.push_back(stat_connectedness);
        result.push_back(stat_hom);
        result.push_back(stat_cohom);
        result.push_back(stat_hom_cohom);
        return result;
    }

    // Get a pair of paths with error delta: sortest path / connectedness path
    std::pair<std::vector<Operation>, std::vector<Operation> > get_pair_paths_delta(const Hdvf_type& X, int delta) {
        std::pair<std::vector<Operation>, std::vector<Operation> > res;
        std::vector<Operation> path, path_connected;
        bool found=false;
        for (Map_type::const_iterator it = map.cbegin(); !found && (it != map.cend()); ++it) {
            if ((it->second.distance_connectedness - it->second.distance) == delta) {
                found = true;
                HDVF_type X1(build_hdvf_from_psc_flags(it->first, dim));
                // Path for connectedness
                path_connected = connectedness(X, X1, dim);
                // Reconstruct flooding path (shortest path)
                size_t id_X(map.at(X.psc_flags(dim)).id), id_current(it->second.id);
                while (id_current != id_X) {
                    path.push_back(map.at(id_to_flags.at(id_current)).op);
                    id_current= map.at(id_to_flags.at(id_current)).pred;
                }
                std::reverse(path.begin(), path.end());
                res.first=path;
                res.second=path_connected;
                return res;
            }
        }
        // If delta not found, raise an error
        std::cerr << "Error : delta not found" << std::endl;
        throw "get_pair_paths_delta error";
    }

    void export_path_to_vtk(const Hdvf_type& X, const std::vector<Operation>& path, std::string root_name) {
        CGAL::IO::write_VTK(hdvf, complex, (root_name+"_0"));
        Hdvf_type X_current(X);
        for (int i=0; i<path.size(); ++i) {
            path.at(i).afficher();
            X_current = operer(X_current, path.at(i));
            CGAL::IO::write_VTK(X_current, complex, (root_name+"_"+to_string(i+1)));
        }
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

    std::vector<size_t> get_min_hom_generators (int min_hom_len) {
        std::vector<size_t> res;
        for (Map_type::const_iterator it = map.cbegin(); it != map.cend(); ++it) {
            if ((it->second.len_hom) == min_hom_len)
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
    Hdvf_space<HDVF_type> hs(complex, 1, filename);
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

//        std::cout << " === three_triangles =============================== " << std::endl;
//        tab.push_back("data/simp/three_triangles_f.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        tab.push_back("data/simp/three_triangles_1_1.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        tab.push_back("data/simp/three_triangles_2_1.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        tab.push_back("data/simp/three_triangles_2_3.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        tab.push_back("data/simp/three_triangles_3_2.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        tab.push_back("data/simp/three_triangles_v.simp");
//        tab_nodes.push_back("data/simp/three_triangles.nodes");
//
//        std::cout << " === six_triangles =============================== " << std::endl;
//
//        tab.push_back("data/simp/six_triangles_f.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_1.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_2.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_3.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_4_1.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_4_2.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_4_3.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_5_1.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/six_triangles_v.simp");
//        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        std::cout << " === larger models =============================== " << std::endl;

//        tab.push_back("data/simp/HDVF_size/three_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/three_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/six_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/six_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/seven_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/seven_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/height_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/height_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/nine_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/nine_triangles.nodes");

//        tab.push_back("data/simp/HDVF_size/ten_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/ten_triangles.nodes");

        std::cout << " === stats delta =============================== " << std::endl;

        //        tab.push_back("data/simp/six_triangles_f.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
                tab.push_back("data/simp/six_triangles_1.simp");
                tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
        //        tab.push_back("data/simp/six_triangles_2.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
        //        tab.push_back("data/simp/six_triangles_3.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
        //        tab.push_back("data/simp/six_triangles_4_1.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
        //        tab.push_back("data/simp/six_triangles_5_1.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
        //
        //        tab.push_back("data/simp/six_triangles_v.simp");
        //        tab_nodes.push_back("data/simp/six_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/three_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/three_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/six_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/six_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/seven_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/seven_triangles.nodes");
//
//        tab.push_back("data/simp/HDVF_size/height_triangles_3.simp");
//        tab_nodes.push_back("data/simp/HDVF_size/height_triangles.nodes");

        for (int i=0; i<tab.size(); ++i) {
            compute_stat(tab.at(i), tab_nodes.at(i));
        }
    }

    return 0;
}
