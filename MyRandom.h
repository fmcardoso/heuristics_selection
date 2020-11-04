//
// Created by fmcardoso on 9/27/19.
//
#include <vector>
#include <deque>

#ifndef COOPGA11_MYRANDOM_H
#define COOPGA11_MYRANDOM_H

#define NormRANu (2.3283063671E-10F)

class MyRandom {
public:
    MyRandom(int seed);
    double Random();
private:
    unsigned  int irr[256];
    unsigned int ir1;
    unsigned char ind_ran,ig1,ig2,ig3;
    void init_ran(int);
};

class uniform_real_distribution{
public:
    uniform_real_distribution(double a, double b);
    double Random(MyRandom *rng);
private:
    double a, b;
};

class uniform_int_distribution{
public:
    uniform_int_distribution(int a, int b);
    int Random(MyRandom *rng);
private:
    int a, b;
};

class uniform_int8_distribution{
public:
    uniform_int8_distribution(int8_t a, int8_t b);
    int8_t Random(MyRandom *rng);
private:
    int8_t a, b;
};

class RandomHelper{
public:

/**
 * Fischer-Yates shuffling as specified by Knuth.
 * @tparam T
 * @tparam A
 * @param vect
 * @param rng
 */
    template<typename T, typename A>
    static void shuffle(std::vector<T, A> *vect, MyRandom *rng) {
        int size = vect->size();
        int aux, j;

        for (int i =0; i < (size - 1); i ++){
            uniform_int_distribution indexDist(i, vect->size() -1);
            aux = (*vect)[i];
            j = indexDist.Random(rng);
            (*vect)[i] = (*vect)[j];
            (*vect)[j] = aux;
        }
    }

    template<typename T, typename A>
    static void shuffle(std::deque<T, A> *vect, MyRandom *rng) {
        int size = vect->size();
        int aux, j;

        for (int i =0; i < (size - 1); i ++){
            uniform_int_distribution indexDist(i, vect->size() -1);
            aux = (*vect)[i];
            j = indexDist.Random(rng);
            (*vect)[i] = (*vect)[j];
            (*vect)[j] = aux;
        }
    }

    };




#endif //COOPGA11_MYRANDOM_H
