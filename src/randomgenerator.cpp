#include "randomgenerator.hpp"

#include <cfloat>  // DBL_MAX
#include <cmath>   // std::nextafter
#include <iostream>
#include <random>

using namespace std;

// Helpers
RandomGenerator::RandomGenerator(long long int s)
    : generator_normal(s), generator_uniform(s), generator_uniform_long(s) {}

double RandomGenerator::randomNormal(double mean, double std) {
    static std::normal_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(generator_normal, parm_t{mean, std});
}

double RandomGenerator::random48() {
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(generator_uniform, parm_t{0.0, std::nextafter(1.0, DBL_MAX)});
}

double RandomGenerator::randomRange(const double from, const double thru) {
    return std::uniform_real_distribution<double>{from,
                                                  thru}(generator_uniform);
}

int RandomGenerator::randomRange(const int from, const int thru) {
    return std::uniform_int_distribution<int>{from, thru}(generator_uniform);
}

long long int RandomGenerator::randomRangeLong(const long long int from,
                                               const long long int thru) {
    return std::uniform_int_distribution<long long int>{
        from, thru}(generator_uniform_long);
}
