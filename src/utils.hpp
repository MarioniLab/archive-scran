#ifndef UTILS_HPP
#define UTILS_HPP
#include <cmath>

struct Rx_random_seed {
    Rx_random_seed() {
        GetRNGstate();
    }
    ~Rx_random_seed() {
        PutRNGstate();
    }
};

template <typename Pt>
void Rx_shuffle (Pt start, Pt end) {
    size_t n=(end-start), i;
    for (i=n-1; i>0; --i) {
        std::swap(start[i], start[int(std::floor(unif_rand()*(i+1)))]); // Equal chance of getting anything from [0, i].
    }
    return;
}

#endif
