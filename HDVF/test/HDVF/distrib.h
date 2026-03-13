//
//  distrib.hpp
//  test1
//
//  Created by Bac Alexandra on 21/03/2022.
//

#ifndef distrib_hpp
#define distrib_hpp

#include <iostream>
using namespace std ;

template <typename T>
class distrib {
    vector<T> _data ;
    T _mean_sum, _min, _max, _stdev ;
    size_t _number_of_bins; // If 0, bins have size 1 (for int data)
public:
    distrib() : _mean_sum(0.), _min(0.), _max(0.), _stdev(0.), _number_of_bins(0)  { } ;
    distrib(const distrib<T>& d, size_t number_of_bins = 0) : _data(d._data), _mean_sum(d._mean_sum), _min(d._min), _max(d._max), _stdev(d._stdev), _number_of_bins(number_of_bins) {}
    distrib(vector<T> &vec) : _data(vec), _mean_sum(0.), _min(0.), _max(0.), _stdev(0.), _number_of_bins(vec._number_of_bins) {
        if (_data.size()>0)
        {
            _mean_sum = _data.at(0) ;
            _min = _data.at(0) ;
            _max = _data.at(0) ;
            for (int i=1; i<_data.size(); ++i)
            {
                _mean_sum += _data.at(i) ;
                if (_data.at(i) < _min)
                    _min = _data.at(i) ;
                if (_data.at(i) > _max)
                    _max = _data.at(i) ;
            }
        }
        update_stdev() ;
    } ;
    
    template <typename TT>
    friend std::ostream& operator<<(std::ostream& out, distrib<TT>& dist);
    
    void update_stdev ()
    {
        T mean = _mean_sum / _data.size(), tmp ;
        for (int i=0; i<_data.size(); ++i)
        {
            tmp = _data.at(i) - mean ;
            _stdev += tmp*tmp ;
        }
        _stdev = sqrt(_stdev / _data.size()) ;
    } ;
    // Accesseurs
    T get_mean () const { return _mean_sum/_data.size() ; } ;
    T get_stdev () { update_stdev() ; return _stdev ; } ;
    T get_min () const { return _min ; } ;
    T get_max () const { return _max ; } ;
    std::vector<size_t> get_hist () {
        std::vector<size_t> res;
        if (_number_of_bins == 0)
            res.resize(_max+1);
        else
            res.resize(_number_of_bins);

        for (T data : _data) {
            // Compute the index of the bin
            size_t i;
            if (_number_of_bins == 0)
                i=data;
            else
                i = T(((data-_min)*_number_of_bins) / (_max-_min));
            ++res.at(i);
        }
        return res;
    }
    std::vector<T> get_hist_labels() {
        std::vector<T> res;
        if (_number_of_bins == 0) {
            for (T i=0; i<=_max; ++i)
                res.push_back(i);
        }
        else {
            res.push_back(_min);

            for (int i=1; i<= _number_of_bins; ++i) {
                T tmp((_min*_number_of_bins + i*(_max-_min))/_number_of_bins);
                res.push_back(tmp);
            }
        }
        return res;
    }

    void add_data (T dat)
    {
        if (_data.size() == 0)
        {
            _min = dat ;
            _max = dat ;
        }
        else
        {
            if (dat < _min)
                _min = dat ;
            if (dat > _max)
                _max = dat ;
        }
        _data.push_back(dat) ; _mean_sum += dat ;
    }
    void clear ()
    { _data.clear() ; _min = 0. ; _max = 0. ; _mean_sum = 0. ; _stdev = 0. ;}
    void to_matlab () {
        std::cout << "A = [" ;
        for (T d : _data)
            std::cout << d << "; ";
        std::cout << "] ;" << std::endl;
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, distrib<T>& dist) {
    out << "K min : " << dist.get_min() << " - max : " << dist.get_max() << std::endl ;
    out << "K mean : " << dist.get_mean() << " - sigma : " << dist.get_stdev() << std::endl;
    return out;
}

#endif /* distrib_hpp */
