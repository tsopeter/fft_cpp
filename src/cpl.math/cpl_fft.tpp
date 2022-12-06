#include "cpl_fft.hpp"

#define FP_PI 3.1415

cpl_fft_implt::IQElement WNKn(size_t N, double n, double k, double c = -1.0) {
    double argument = 2.0 * FP_PI * n * k / static_cast<double>(N);
    double real = std::cos(argument);
    double imag = std::sin(argument);
    return {real, c * imag};
}

std::pair<double, double> cpl_fft_implt::IQElement::string_parser (const std::string complex) {
    // We separate into two regions "real"+/-j"imag"
    std::string sreal, simag;

    double negreal = 1;
    double negimag = 1;
    bool left = true;

    for (const char &c : complex) {
        if (c != '+' && c != '-') {
            if (c == 'j') left = false;
            else
                if (left)
                    sreal += c;
                else
                    simag += c;
        }
        else if (c == '-'){
            if (left)
                negreal = -1;
            else
                negimag = -1;
        }
    }

    if (sreal.length() <= 0) sreal = "0.0";
    if (simag.length() <= 0) simag = "0.0";

    std::pair<double, double> ret = {
        std::stod(sreal) * negreal,
        std::stod(simag) * negimag
    };

    return ret;
}

template <size_t N>
cpl_fft_implt::FFT<N>::FFT () {
    // populate the fft_matrix
    for (size_t k = 0; k < N; ++k) {
        for (size_t n = 0; n < N; ++n) {
            IQElement e = WNKn (N, n, k, -1.0);
            IQElement b = WNKn (N, n, k, 1.0);
            fft_matrix[k][n] = e;
            ifft_matrix[k][n] = b;
        }
    }
}

template <size_t N>
cpl_fft_implt::FFT<N>::~FFT () {

}

template <size_t N>
void cpl_fft_implt::FFT<N>::set_input (const std::array<IQElement, N> &_input) {
    for (size_t i  = 0; i < N; ++i)
        input[i] = _input[i];
}

template <size_t N>
size_t cpl_fft_implt::FFT<N>::compute_fft () {
    // multiply the input matrix by fft_matrix
    size_t mcount = 0;
    for (size_t k = 0; k < N; ++k) {
        IQElement acc {};
        for (size_t n = 0; n < N; ++n) {
            const IQElement &b = fft_matrix[k][n];
            const IQElement &x = input[n];

            // do complex multiplication and sum
            acc = acc + (b * x);
        }
        mcount += acc.mcount;
        output[k] = acc;
    }
    return mcount;
}

template <size_t N>
size_t cpl_fft_implt::FFT<N>::compute_ifft () {
    size_t mcount = 0;
    for (size_t n = 0; n < N; ++n) {
        IQElement acc {};
        for (size_t k = 0; k < N; ++k) {
            const IQElement &X = output[k];
            const IQElement &b = ifft_matrix[k][n];
            
            acc = acc + (X * b);
        }
        mcount += acc.mcount;
        output[n] = acc;
    }
    return mcount;
}