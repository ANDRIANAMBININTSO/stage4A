#ifndef CGAL_HDVF_HDVF_SPACE_H
#define CGAL_HDVF_HDVF_SPACE_H


#include <CGAL/license/HDVF.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <random>
#include <CGAL/OSM/Bitboard.h>
#include <CGAL/HDVF/Hdvf_core.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>


namespace CGAL {
namespace Homological_discrete_vector_field {

template<typename ChainComplex>
class Hdvf_space : public Hdvf<ChainComplex> {
    /*! \brief Chain complex type */
    typedef ChainComplex Chain_complex;

    /*! \brief Type of coefficients used to compute homology. */
    typedef typename Chain_complex::Coefficient_ring Coefficient_ring;


    /*!
     Type of parent Hdvf_core class.
     */
    typedef Hdvf_core<Chain_complex, OSM::Sparse_chain, OSM::Sparse_matrix> Base ;
    //typedef Hdvf<Chain_complex> Base;

    // Inherited types
    using Column_chain = typename Base::Column_chain;
    using Row_chain = typename Base::Row_chain;
    using Column_matrix = typename Base::Column_matrix;
    using Row_matrix = typename Base::Row_matrix;


public:
    using Hdvf<Chain_complex>::hd;
    //using Hdvf<Chain_complex>::htdt;

    Hdvf_space(const Chain_complex& K, int hdvf_opt = OPT_FULL, int dimension_restriction = -1): Hdvf<Chain_complex>(K){ }
    //Hdvf_space(const Hdvf_space& hdvf) : Hdvf_core<ChainComplex, OSM::Sparse_chain, OSM::Sparse_matrix>(hdvf) { }
    Hdvf_space(const Hdvf_space& hdvf) : Hdvf<Chain_complex>(hdvf) {}
    ~Hdvf_space() { }
};
}
}
#endif //CGAL_HDVF_HDVF_SPACE_H