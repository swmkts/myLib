#ifndef DIRICHLET_HPP__
#define DIRICHLET_HPP__

#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <random>
#include <cassert>

namespace myLib {
    namespace stat {

template <class RealType = double, class VecType = std::vector<double> >
class dirichlet_distribution
{
public:
    typedef RealType value_type;
    typedef VecType result_type;

    dirichlet_distribution(VecType alpha = VecType{1.0, 1.0, 1.0})
        : _alpha(alpha)
    {

    }

    VecType alpha() const
    {
        return _alpha;
    }

private:
    VecType _alpha;

};

typedef dirichlet_distribution<double> dirichlet;

template <class RealType = double, class VecType = std::vector<RealType> >
inline VecType mean(const dirichlet_distribution<RealType>& dist)
{
    auto mean = dist.alpha();
    RealType sum = std::accumulate(mean.begin(), mean.end(), 0.0);
    std::for_each(mean.begin(), mean.end(), [sum](RealType& v) { v /= sum; });

    return mean;
}

template <class RealType = double, class VecType = std::vector<RealType> >
inline VecType var(const dirichlet_distribution<RealType>& dist)
{
    auto var = dist.alpha();
    RealType sum = std::accumulate(var.begin(), var.end(), 0.0);
    std::for_each(var.begin(), var.end(), [sum](RealType& v) { v = v * (sum - v) / (sum * sum * (sum + 1)); });

    return var;
}


template <class RealType, class VecType = std::vector<RealType> >
inline RealType pdf(const dirichlet_distribution<RealType>& dist, const VecType& x)
{
    auto alpha = dist.alpha();
    RealType sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
    auto gamma_alpha = alpha;
    std::for_each(gamma_alpha.begin(), gamma_alpha.end(), tgamma);
    RealType dens = tgamma(sum) / std::accumulate(gamma_alpha.begin(), gamma_alpha.end(), 1.0, std::multiplies<RealType>());
    for (int i = 0, max_i = alpha.size(); i < max_i; ++i) {
        dens *= pow(x[i], alpha[i] - 1);
    }
    return dens;
}

    } // myLib::stat

} // myLib

#endif // DIRICHLET_HPP__

