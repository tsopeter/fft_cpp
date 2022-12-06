#include <iostream>
#include "cpl.math/cpl.math.hpp"

template <size_t N>
bool validate (const cpl_math::cint errate, const std::array<cpl_math::cint, N> &a, const std::array<cpl_math::cint, N> &b) {
    bool valid = true;
    for (size_t i = 0; i < N; ++i) {
        const cpl_math::cint &A = a[i];
        const cpl_math::cint &B = b[i];

        const cpl_math::cint &C = A - B;
        if (C > errate)
            valid = false;
    }
    return valid;
}

template <size_t N = 9>
std::array<cpl_math::cint, N> nine_point(std::array<cpl_math::cint, N> xs) {

    std::array<std::pair<cpl_math::cint, cpl_math::cint>, 9> NN9;
    std::array<std::pair<cpl_math::cint, cpl_math::cint>, 3> NN3;

    for (size_t i = 0; i < N; ++i) {
        NN9[i].first = WNKn(N, 1, i);
        NN9[i].second = WNKn(N, 2, i);
    }

    for (size_t i = 0; i < N/3; ++i) {
        NN3[i].first = WNKn(N/3, 1, i);
        NN3[i].second = WNKn(N/3, 2, i);
    }

    std::array<cpl_math::cint, 3> A = {
        xs[0] + xs[3] * NN3[0].first + xs[6] * NN3[0].second,
        xs[0] + xs[3] * NN3[1].first + xs[6] * NN3[1].second,
        xs[0] + xs[3] * NN3[2].first + xs[6] * NN3[2].second
    };

    std::array<cpl_math::cint, 3> B = {
        xs[1] + xs[4] * NN3[0].first + xs[7] * NN3[0].second,
        xs[1] + xs[4] * NN3[1].first + xs[7] * NN3[1].second,
        xs[1] + xs[4] * NN3[2].first + xs[7] * NN3[2].second
    };

    std::array<cpl_math::cint, 3> C = {
        xs[2] + xs[5] * NN3[0].first + xs[8] * NN3[0].second,
        xs[2] + xs[5] * NN3[1].first + xs[8] * NN3[1].second,
        xs[2] + xs[5] * NN3[2].first + xs[8] * NN3[2].second
    };

    std::array<cpl_math::cint, 9> Xs = {
        A[0] + B[0] * NN9[0].first + C[0] * NN9[0].second,
        A[1] + B[1] * NN9[1].first + C[1] * NN9[1].second,
        A[2] + B[2] * NN9[2].first + C[2] * NN9[2].second,
        A[0] + B[0] * NN9[3].first + C[0] * NN9[3].second,
        A[1] + B[1] * NN9[4].first + C[1] * NN9[4].second,
        A[2] + B[2] * NN9[5].first + C[2] * NN9[5].second,
        A[0] + B[0] * NN9[6].first + C[0] * NN9[6].second,
        A[1] + B[1] * NN9[7].first + C[1] * NN9[7].second,
        A[2] + B[2] * NN9[8].first + C[2] * NN9[8].second
    };

    return Xs;
}

auto main()->int{
    cpl_math::cint e = {0.5, 3.2};
    e = e^(-3);

    std::cout<<e<<'\n';

    constexpr size_t N = 9;
    cpl_math::fft_core<N> fft {};

    const std::array<cpl_math::cint, N> xs = {
        "1.0+j1.5",
        "2.0-j0.5",
        "3.0+j3.2",
        "-4.0-j5.5",
        "5.0+j6.1",
        "6.0+j7.1",
        "-7.0-j5.56",
        "8.0+j10.4",
        "9.0-j6.6"
    };

    fft.set_input(xs);
    auto x = fft.compute_fft();
    std::cout<<fft<<'\n';

    auto rr = nine_point(xs);

    std::cout<<"It took matrix FFT: "<<x<<" multiply operations\n";
    for (const auto &r : rr) {
        std::cout<<r<<'\n';
    }

    std::cout<<(validate<9>({0.01, 0.01}, rr, fft.get_output()) ? "Valid!" : "Invalid!")<<'\n';

    return 0;
}