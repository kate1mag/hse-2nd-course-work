/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testadvcd.cpp
 * Author: posypkin
 *
 * Created on June 28, 2019, 5:03 PM
 */
#include <iostream>
#include <iterator>
#include "advancedcoordescent.hpp"

constexpr int n = 3;

//double f(const double* x) {
//    double v;
//    for (int i = 0; i < n; i++)
//        v += x[i] * x[i];
//    return v;
//}

struct F {

    double operator()(const double *x) {
        double v;
        mCnt ++;
        for (int i = 0; i < n; i++)
            v += x[i] * x[i];
        return v;
    }
    
    int mCnt = 0;
};

/*
 * 
 */
int main(int argc, char** argv) {

    panther::AdvancedCoorDescent<double> adv;
    double x[n];
    double a[n], b[n];
    std::fill(a, a + n, -1);
    std::fill(b, b + n, 2);
    std::fill(x, x + n, 1);
    F f;
    double v = adv.search(n, x, a, b, std::ref(f));
    std::cout << "Found " << v << " at [";
    std::copy(x, x + n, std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    std::cout << f.mCnt << " function calls done\n";


    return 0;
}

