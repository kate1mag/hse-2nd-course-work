/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <iterator>
#include "gridlip.hpp"

constexpr int n = 3;

double f(const double* x) {
    double v;
    for(int i = 0; i < n; i ++)
        v += x[i] * x[i];
    return v;
}

int main() {
    panther::GridLip<double> gridlip;
    gridlip.mOptions.mEps = 1e-3;
    double x[n];
    double a[n], b[n];
    std::fill(a, a + n, -1.01);
    std::fill(b, b + n, 2.57);
    double v = gridlip.search(n, x, a, b, f);
    std::cout << "Found " << v << " at [" ;
    std::copy(x, x + n, std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
}