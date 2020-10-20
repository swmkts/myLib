#ifndef DIRICHLET_DISTRIBUTION_HPP__
#define DIRICHLET_DISTRIBUTION_HPP__

#include <iostream>
#include <istream>
#include <ostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <cassert>

namespace myLib {
    namespace stat {

/** dirichlet distribution ramdom variable generator class
 * boost-like implementation.
 **/
template<class RealType = double, class VecType= std::vector<RealType> >
class dirichlet_distribution
{
public:
    typedef RealType value_type;
    typedef VecType result_type;

    class param_type
    {
    public:
        typedef dirichlet_distribution distribution_type;

        param_type(VecType& alpha_arg = VecType{1.0, 1.0})
            : _alpha(alpha_arg) {}
        
        VecType alpha() const
        {
            return _alpha;
        }

        
    
    private:
        VecType _alpha;
    };

    dirichlet_distribution(const VecType& alpha_arg = VecType{1.0, 1.0, 1.0})
        : _alpha(alpha_arg)
    {

    }

    VecType alpha() const { return _alpha; }

    param_type param() const { return param_type(_alpha); }

    void param(const param_type& parm)
    {
        _alpha = parm.alpha();
    }

    void reset() { }

    /**
     * Returns dirichlet distributed vector.
     * By normalization of gamma distributed vector.
     **/
    template<class Engine>
    result_type operator()(Engine& eng)
    {
        std::vector<double> dirichlet_rv;
        double sum_dir = 0;
        for (int i = 0, alpha_len = _alpha.size(); i < alpha_len; ++i) {
            std::gamma_distribution<double> g_i(_alpha[i], 1.0);
            double dir_i = g_i(eng);
            dirichlet_rv.push_back(dir_i);
            sum_dir += dir_i;
            std::cout << dir_i << std::endl;
        }

        std::for_each(dirichlet_rv.begin(), dirichlet_rv.end(), [sum_dir](double& v) { v /= sum_dir; });

        return dirichlet_rv;
    }

    template<typename CharT, typename Traits>
    friend std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
                const dirichlet_distribution<RealType>& x)
    {
        std::vector<RealType> alpha = x.alpha();
        for (auto al: alpha) {
            os << al << " ";
        }
        return os;
    }

    // TODO: vector alpha input operator
    // template<typename CharT, typename Traits>
    // friend std::basic_istream<CharT, Traits>&
    // opeator>>(std::basic_istream<CharT, Traits>& is,
    //             const dirichlet_distribution<RealType>& x)
    // {
        
    // }


private:
    VecType _alpha;
};

} // myLib::stat

} // myLib

#endif // DIRICHLET_DISTRIBUTION_HPP__

