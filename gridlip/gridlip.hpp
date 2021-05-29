/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   gridlip.hpp
 * Author: posypkin
 *
 * Created on July 13, 2019, 1:28 PM
 */

#ifndef GRIDLIP_HPP
#define GRIDLIP_HPP


#include <math.h>
#include <algorithm>
#include <vector>
#include <limits>
#include <common/bbsolver.hpp>

/**
 * Simple grid-based Lipschitzian solver developed by 
 * Alexey Gorbunov and Mikhail Posypkin
 */

namespace panther {

    /* hyperinterval info */
    template <class T>
    class Box {
        int size;
    public:
        T *mA = nullptr;
        T *mB = nullptr;
        T mLocUB, mLocLO;

        Box(int n, const T* a, const T* b) {
            try {
                size = n;
                this->mA = new T[size];
                this->mB = new T[size];
            } catch (std::bad_alloc e) {
                throw e;
            }
            for (int i = 0; i < size; i++) {
                this->mA[i] = a[i];
                this->mB[i] = b[i];
            }
        }

        Box(Box&& p) {
            std::swap(mA, p.mA);
            std::swap(mB, p.mB);
            this->mLocLO = p.mLocLO;
            this->mLocUB = p.mLocUB;
        }

        Box & operator=(Box && p) {
            std::swap(mA, p.mA);
            std::swap(mB, p.mB);
            this->mLocLO = p.mLocLO;
            this->mLocUB = p.mLocUB;
            return *this;
        }

        void change(const T*a, const T*b) {
            for (int i = 0; i < this->size; i++) {
                this->mA[i] = a[i];
                this->mB[i] = b[i];
            }
        }

        ~Box() {
            delete [] mA;
            delete [] mB;
        }
    };

    /**
     * A simple black-box optimizer that uses Lipschitzian bounds
     */
    template <class T> class GridLip : public BlackBoxSolver <T> {
    public:

        struct Options {
            // Accuracy
            T mEps = 1e-1;
            // Number of node per dimension
            int mNodes = 4;
        } mOptions;

        /**
         * Constructor
         */
        GridLip() {
        }

        ~GridLip() {
            if (step != nullptr) delete[]step;
            if (x != nullptr) delete[]x;
            if (Fvalues != nullptr) delete[]Fvalues;
        }

        /**
         * Search with grid solver
         * @param n number of task dimensions
         * @param x coordinates of founded minimum (retvalue)
         * @param a,b left/right bounds of search region
         * @param f pointer to function for which search minimum
         */
        virtual T search(int n, T* xfound, const T * const a, const T * const b, const std::function<T(const T * const)> &f) {
            /* reset variables */
            dim = n;
            nodes = mOptions.mNodes;
            eps = mOptions.mEps;
            allnodes = static_cast<int> (pow(nodes, dim));
            /* create 2 vectors */
            /* P contains parts (hyperintervals on which search must be performed */
            /* P1 temporary */
            try {
                /* step of grid in every dimension */
                if (step != nullptr) delete[]step;
                if (x != nullptr) delete[]x;
                if (Fvalues != nullptr) delete[]Fvalues;
                step = new T[dim];
                x = new T[dim];
                Fvalues = new T[allnodes];
            } catch (std::bad_alloc& ba) {
                std::cerr << ba.what() << std::endl;
                return UPB;
            }
            std::vector<Box <T> > P, P1;
            /* Upper bound */
            UPB = std::numeric_limits<T>::max();

            /* Add first hyperinterval */

            try {
                P.emplace_back(dim, a, b);
            } catch (std::exception& e) {
                std::cerr << e.what() << std::endl;
                return UPB;
            }

            /* If hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates */
            /* xs temporary array for storage local min coordinates */
            T *a1, *b1, *xs;
            try {
                a1 = new T[dim];
                b1 = new T[dim];
                xs = new T[dim];
            } catch (std::bad_alloc& ba) {
                std::cerr << ba.what() << std::endl;
                return UPB;
            }

            /* Each hyperinterval can be subdivided or pruned (if non-promisable or fits accuracy) */
            while (!P.empty()) {
                /* number of iterations on this step (BFS) */
                unsigned int parts = P.size();

                /* For all hyperintervals on this step perform grid search */
                for (unsigned int i = 0; i < parts; i++) {
                    /* local values of upper and lower bounds, value of delta*L (Lipshitz const) */
                    T lUPB, lLOB, ldeltaL;
                    T* ta = P[i].mA, *tb = P[i].mB;
                    gridEvaluator(ta, tb, xs, &lUPB, &lLOB, &ldeltaL, f);
                    P[i].mLocLO = lLOB;
                    P[i].mLocUB = lUPB;
                    /* remember new results if less then previous */
                    updateRecords(lUPB, xfound, xs);
                }

                /* Choose which hyperintervals should be subdivided */
                for (unsigned int i = 0; i < parts; i++) {
                    /* Subdivision criteria */
                    if (P[i].mLocLO < (UPB - eps)) {
                        /* If subdivide, choose dimension (the longest side) */
                        int choosen = chooseDim(P[i].mA, P[i].mB);

                        /* Make new edges for 2 new hyperintervals */
                        for (int j = 0; j < dim; j++) {
                            if (j != choosen) { /* [a .. b1] [a1 .. b] */
                                a1[j] = P[i].mA[j]; /* where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ] */
                                b1[j] = P[i].mB[j]; /* and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ] */
                            } else {
                                a1[j] = P[i].mA[j] + fabs(P[i].mB[j] - P[i].mA[j]) / 2.0;
                                b1[j] = a1[j];
                            }
                        }
                        /* Add 2 new hyperintervals, parent HI no longer considered */
                        try {
                            P1.emplace_back(dim, P[i].mA, b1);
                        } catch (std::exception& e) {
                            std::cerr << e.what() << std::endl;
                            return UPB;
                        }
                        try {
                            P1.emplace_back(dim, a1, P[i].mB);
                        } catch (std::exception& e) {
                            std::cerr << e.what() << std::endl;
                            return UPB;
                        }
                    }
                }

                P.clear();
                P.swap(P1);
            }
            delete[]a1;
            delete[]b1;
            delete[]xs;
            P.clear();
            return UPB;
        }


    private:

        T eps; /* required accuracy */
        int nodes, dim, allnodes; /* internal varibale for handlig errors and number of nodes per dimension */
        T UPB, LOB; /* obtained upper bound and lower bound */
        T *x = nullptr, *step = nullptr, *Fvalues = nullptr;

        /* Get R (reliable coefficient) for the corresponding step lenght*/
        virtual double getR(const T delta) {
            /* The value depends on step lenght (test formula, may be changed) */
            return static_cast<double> (exp(delta));
        }

        /* Update the current record and its coordinates in accordance with new results obtained on some hyperinterval*/
        virtual void updateRecords(const T LU, T* x, const T *xs) {
            if (LU < UPB) {
                UPB = LU;
                for (int i = 0; i < dim; i++) {
                    x[i] = xs[i];
                }
            }
        }

        /* Select dimension for subdivide hyperinteral (choose the longest side) */
        virtual int chooseDim(const T *a, const T *b) {
            T max = std::numeric_limits<T>::min(), cr;
            int i, maxI = 0;
            for (i = 0; i < dim; i++) {
                cr = fabs(b[i] - a[i]);
                if (cr > max) {
                    max = cr;
                    maxI = i;
                }
            }
            return maxI;
        }

        virtual void gridEvaluator(const T *a, const T *b, T* xfound, T *Frp, T *LBp, T *dL, const std::function<T(const T * const)> &compute) {
            T Fr = std::numeric_limits<T>::max(), L = std::numeric_limits<T>::min(), delta = 0, LB;
            double R;
            for (int i = 0; i < dim; i++) {
                step[i] = fabs(b[i] - a[i]) / (nodes - 1);
                delta += 0.5 * step[i];
            }
            R = getR(delta);
            int node = 0;
            /* Calculate and cache the value of the function in all points of the grid */
            for (int j = 0; j < allnodes; j++) {
                int point = j;
                for (int k = dim - 1; k >= 0; k--) {
                    int t = point % nodes;
                    point = (int) (point / nodes);
                    x[k] = a[k] + t * step[k];
                }
                T rs = compute(x);
                Fvalues[j] = rs;
                /* also remember minimum value across the grid */
                if (rs < Fr) {
                    Fr = rs;
                    node = j;
                }
            }


            /* Calculate coordinates of obtained upper bound */

            for (int k = dim - 1; k >= 0; k--) {
                int t = node % nodes;
                node = (int) (node / nodes);
                xfound[k] = a[k] + t * step[k];
            }

            /* Calculate all estimations of Lipshitz constant and choose maximum estimation */
            for (int j = 0; j < allnodes; j++) {
                int neighbour;
                for (int k = 0; k < dim; k++) {
                    int board = (int) pow(nodes, k + 1);
                    neighbour = j + board / nodes;
                    if ((neighbour < allnodes) && ((j / board) == (neighbour / board))) {
                        T loc = fabs(Fvalues[j] - Fvalues[neighbour]) / step[dim - 1 - k];
                        L = loc > L ? loc : L;
                    }
                }
            }

            /* final calculation */

            LB = R * L * delta;
            *dL = LB;
            LB = Fr - LB;
            *Frp = Fr;
            *LBp = LB;
        }
    };
}

#endif /* GRIDLIP_HPP */

