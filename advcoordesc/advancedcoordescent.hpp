/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   advancedcoordescent.hpp
 * Author: posypkin
 *
 * Created on June 28, 2019, 4:59 PM
 */

#ifndef ADVANCEDCOORDESCENT_HPP
#define ADVANCEDCOORDESCENT_HPP

#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>
#include <common/bbsolver.hpp>
#include <common/vec.hpp>

namespace panther {

    /**
     * An adaptive advanced coordinate descent solver
     */
    template <class T> class AdvancedCoorDescent : public BlackBoxSolver <T> {
    public:

        struct Options {
            // Initial step (positive)
            T mInitStep = 1e-1;
            // Incremental parameter
            T mInc = 1.5;
            // Decremental parameter
            T mDec = 0.5;
            // Minimal step size
            T mMinStep = 1e-3;
        } mOptions;

        T search(int n, T* x, const T * const a, const T * const b, const std::function<T(const T * const)> &f) override {
            std::vector<T> sft(n, mOptions.mInitStep);
            auto maxStep = [&sft, n]() {
                T rv = 0;
                for (int i = 0; i < n; i++) {
                    rv = std::max(rv, sft[i]);
                }
                return rv;
            };
            T v = f(x);
            while (maxStep() >= mOptions.mMinStep) {
                for (int i = 0; i < n; i++) {
                    const T h = sft[i];
                    const T xi = x[i];
                    x[i] = std::min(x[i] + h, b[i]);
                    T vn = f(x);
                    if (vn < v) {
                        v = vn;
                        sft[i] *= mOptions.mInc;
                    } else {
                        x[i] = std::max(x[i] - 2 * h, a[i]);
                        vn = f(x);
                        if (vn < v) {
                            v = vn;
                            sft[i] *= mOptions.mInc;
                        } else {
                            x[i] = xi;
                            sft[i] *= mOptions.mDec;
                        }
                    }
                }
            }
            return v;
        }

    private:


    };
}



#endif /* ADVANCEDCOORDESCENT_HPP */

