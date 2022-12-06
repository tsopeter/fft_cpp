#ifndef cpl_fft_hpp
#define cpl_fft_hpp

#include <iostream>
#include <unistd.h>
#include <string>
#include <array>
#include <cmath>

namespace cpl_fft_implt {
    struct IQElement {
        double real;
        double imag;
        size_t mcount = 0;
        size_t acount = 0;

        friend std::ostream &operator<<(std::ostream &os, IQElement element) {
            return os<<element.to_string();
        }

        IQElement(double real = 0.0f, double imag = 0.0f) : real(real), imag(imag) {}
        
        IQElement(const std::string complex) {
            this->operator=(complex);
        }

        IQElement(const char str[]) {
            this->operator=(std::string(str));
        }

        IQElement (const IQElement &other) {
            this->operator=(other);
        }

        void operator=(const std::string complex) {
            auto pair = string_parser(complex);
            real = pair.first;
            imag = pair.second;
        }

        void operator=(const IQElement &other) {
            real = other.real;
            imag = other.imag;
            mcount = other.mcount;
            acount = other.acount;
        }

        bool operator< (IQElement other) const {
            auto a = abs(), b = other.abs();
            return a < b;
        }

        bool operator<= (IQElement other) const {
            auto a = abs(), b = other.abs();
            return a <= b;
        }

        bool operator> (IQElement other) const {
            auto a = abs(), b = other.abs();
            return a > b;
        }

        bool operator>= (IQElement other) const {
            auto a = abs(), b = other.abs();
            return a >= b;
        }

        double abs() const {
            auto s = real * real + imag * imag;
            return std::sqrt(s);
        }

        IQElement operator+ (IQElement other) const {
            IQElement e = {real + other.real, imag + other.imag};
            e.acount = 2 + other.acount + acount;
            e.mcount = other.mcount + mcount;
            return e;
        }

        IQElement operator* (IQElement other) const {
            auto x = real * other.real - imag * other.imag;
            auto y = real * other.imag + imag * other.real;
            IQElement e = {x, y};
            e.mcount = 4 + other.mcount + mcount;
            e.acount = 2 + other.acount + acount;
            return e;
        }

        IQElement operator* (double other) const {
            IQElement e = {real * other, imag * other};
            e.mcount = 2 + mcount;
            e.acount = acount;
            return e;
        }

        IQElement operator- (IQElement other) const {
            IQElement e = {real - other.real, imag - other.imag};
            e.acount = 2 + other.acount + acount;
            e.mcount = other.mcount + mcount;
            return e;
        }

        IQElement operator^ (int power) const {
            if (power == 0) return {1.0f, 0.0f};

            IQElement acc = {1, 0};
            IQElement ppp = {real, imag};
            bool inverse = power < 0;

            power = std::abs(power);
            for (int i = 0; i < power; ++i)
                acc = acc * ppp;

            if (inverse)
                acc = acc.inverse();

            return acc;
        }

        IQElement inverse() const {
            const IQElement self = {real, imag};
            const IQElement comp = {real, -imag};
            const IQElement num = comp;
            const IQElement den = self * comp;
            const double scalar = 1/den.real;
            return num * scalar;
        }

        std::string to_string() const {
            const std::string sreal = std::to_string(real);
            const std::string simag = std::to_string(std::abs(imag));
            const std::string sep  = (imag < 0) ? "-j" : "+j";
            const std::string ret = sreal + sep + simag;
            return ret;
        }

        std::pair<double, double> string_parser(const std::string complex);
    };

    template <size_t N>
    class FFT {
    public:
        FFT();
        ~FFT();

        void set_input(const std::array<IQElement, N> &_input);

        friend std::ostream &operator<<(std::ostream &os, FFT<N> &fft) {
            os<<"input:\n";
            for (size_t i = 0; i < N; ++i)
                os<<"("<<i<<"): "<<fft.input[i]<<'\n';
            os<<"\noutput:\n";
            for (size_t i = 0; i < N; ++i)
                os<<"("<<i<<"): "<<fft.output[i]<<'\n';
            return os;
        }

        size_t compute_fft();
        size_t compute_ifft();

        std::array<IQElement, N> get_output () {
            return output;
        }

    private:
        std::array<IQElement, N> input  = {};
        std::array<IQElement, N> output = {};
        std::array<std::array<IQElement, N>, N> fft_matrix;
        std::array<std::array<IQElement, N>, N> ifft_matrix;
    };
}

#endif