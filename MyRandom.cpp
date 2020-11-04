#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "MyRandom.h"

double Random(void);
void ini_ran(int SEMILLA);

MyRandom::MyRandom(int seed) {
    init_ran(seed);
}

void MyRandom::init_ran(int SEMILLA){
    int INI,FACTOR,SUM,i;

    srand(SEMILLA);

    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

double MyRandom::Random() {
    double r;

    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;

    return r;
}


uniform_real_distribution::uniform_real_distribution(double a, double b) {
//    __glibcxx_assert(a < b);
    this->a = a;
    this->b=b;
}

double uniform_real_distribution::Random(MyRandom *rng) {
    double r = rng->Random();
    return(a + (b-a)*r);
}

int uniform_int_distribution::Random(MyRandom *rng){
    double r = rng->Random();
    int range = b - a + 1;

    int scaled = (r * range) + a;

    return(scaled);
}

uniform_int_distribution::uniform_int_distribution(int a, int b) {
    this->b =b;
    this->a = a;
}

uniform_int8_distribution::uniform_int8_distribution(int8_t a, int8_t b) {
    this->b =b;
    this->a = a;
}

int8_t uniform_int8_distribution::Random(MyRandom *rng) {
    double r = rng->Random();
    int8_t range = b - a + 1;

    int8_t scaled = (r * range) + a;

    return(scaled);
}
