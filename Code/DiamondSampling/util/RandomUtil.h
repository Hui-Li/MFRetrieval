#ifndef RANDOMUTIL_H
#define RANDOMUTIL_H

#include <random>

template<typename engine_t=std::mt19937>
struct random_number_generator : public engine_t { // {{{
    typedef typename engine_t::result_type result_type;

    random_number_generator(unsigned seed=0): engine_t(seed){ }

    result_type randrange(result_type end=engine_t::max()) { return engine_t::operator()() % end; }
    template<class T=double, class T2=double> T uniform(T start=0.0, T2 end=1.0) {
        return std::uniform_real_distribution<T>(start, (T)end)(*this);
    }
    template<class T=double> T normal(T mean=0.0, T stddev=1.0) {
        return std::normal_distribution<T>(mean, stddev)(*this);
    }
    template<class T=int, class T2=T> T randint(T start=0, T2 end=std::numeric_limits<T>::max()) {
        return std::uniform_int_distribution<T>(start, end)(*this);
    }
    template<class RandIter> void shuffle(RandIter first, RandIter last) {
        std::shuffle(first, last, *this);
    }
};
#endif //RANDOMUTIL_H
