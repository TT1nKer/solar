#ifndef SOLAR_NUMERICS_DETAIL_DOPRI5_DENSE_H
#define SOLAR_NUMERICS_DETAIL_DOPRI5_DENSE_H

#ifndef SOLAR_NUMERICS_DOPRI5_H
#error "include solar/numerics/dopri5.h instead"
#endif

#include <array>
#include <cstddef>

namespace solar::numerics::detail {

template <std::size_t N>
std::array<StateN<N>, 4> dopri5_dense_coefficients(
    const std::array<StateN<N>, 7>& stages) {
    constexpr double extension[7][4] = {
        {1.0,
         -8048581381.0 / 2820520608.0,
         8663915743.0 / 2820520608.0,
         -12715105075.0 / 11282082432.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0,
         131558114200.0 / 32700410799.0,
         -68118460800.0 / 10900136933.0,
         87487479700.0 / 32700410799.0},
        {0.0,
         -1754552775.0 / 470086768.0,
         14199869525.0 / 1410260304.0,
         -10690763975.0 / 1880347072.0},
        {0.0,
         127303824393.0 / 49829197408.0,
         -318862633887.0 / 49829197408.0,
         701980252875.0 / 199316789632.0},
        {0.0,
         -282668133.0 / 205662961.0,
         2019193451.0 / 616988883.0,
         -1453857185.0 / 822651844.0},
        {0.0,
         40617522.0 / 29380423.0,
         -110615467.0 / 29380423.0,
         69997945.0 / 29380423.0},
    };

    std::array<StateN<N>, 4> coefficients{};
    for (std::size_t component = 0; component < N; ++component) {
        for (std::size_t power = 0; power < 4; ++power) {
            for (std::size_t stage = 0; stage < 7; ++stage) {
                coefficients[power][component] +=
                    stages[stage][component] *
                    extension[stage][power];
            }
        }
    }
    return coefficients;
}

} // namespace solar::numerics::detail

#endif // SOLAR_NUMERICS_DETAIL_DOPRI5_DENSE_H
