#ifndef cpl_math_hpp
#define cpl_math_hpp

#include "cpl_fft.tpp"

namespace cpl_math {
    typedef cpl_fft_implt::IQElement cint;

    template <size_t N>
    using fft_core = cpl_fft_implt::FFT<N>;
}

#endif