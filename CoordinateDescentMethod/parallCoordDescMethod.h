#include <cmath>
#include <common/bbsolver.hpp>
#include <common/vec.hpp>
#include <functional>
#include <future>


#include <sstream>
#include <vector>
#include <memory>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <thread>
#include <utility>
#include <mutex>
#include <future>

#ifndef PANTHER_MASTER_PARALLCOORDDESCMETHOD_H
#define PANTHER_MASTER_PARALLCOORDDESCMETHOD_H
std::mutex mtx;
namespace panther {
    template <typename T> class parallCoordinateDescentMethod : public BlackBoxSolver<T> {
    public:
        struct Options {
            T delta = 1e-1;
            T mMinStep = 1e-3;
            T mDec = 0.5;
            T mInc = 1.5;
        } mOptions;

        void search_one_direction(std::vector<T>& result_direction, std::vector<T>& index, T start, int direction,
                const std::function<T(const T * const)> &f, T left_border, T right_border, T* starting_point, int n, T h) {
            T result_for_axis[n];

            for (int j = 0; j < n; ++j) {
                result_for_axis[j] = starting_point[j];
            }

            result_for_axis[direction] = std::min(result_for_axis[direction] + h, right_border);
            result_direction[direction] = f(result_for_axis);
            if (result_direction[direction] < start) {
                index[direction] = direction + 1;
                return;
            } else {
                result_for_axis[direction] = std::max(result_for_axis[direction] - 2 * h, left_border);
                result_direction[direction] = f(result_for_axis);
                index[direction] = -(direction + 1);
                return;
            }
        }

        //покоординатный спуск
        T search(int n, T* starting_point, const T * const left_border, const T * const right_border,
                const std::function<T(const T * const)> &f) override {
            std::vector<T> delta_array(n, mOptions.mMinStep);
            auto maxStep = [&delta_array, n]() {
                T rv = 0;
                for (int i = 0; i < n; i++) {
                    rv = std::max(rv, delta_array[i]);
                }
                return rv;
            };
            std::vector<T> result_direction(n);
            std::vector<T> index(n);
            std::vector<std::thread> threads;
            T start = f(starting_point);
            while (maxStep() >= mOptions.mMinStep) {
                for (int i = 0; i < n; ++i) {    // идём по всем направления
                    std::thread t(std::thread(&parallCoordinateDescentMethod::search_one_direction, this,
                            std::ref(result_direction),std::ref(index), start, i, f, left_border[i],
                            right_border[i], starting_point, n, delta_array[i]));
                    threads.push_back(std::move(t));
                }
                for(std::thread &t: threads){
                    if (t.joinable()) {
                        t.join();
                    }
                }
                T min_result = result_direction[0];
                int new_index = 0;
                for (int i = 0; i < n; ++i) {
                    min_result = std::min(result_direction[i], min_result);
                    delta_array[i] = result_direction[i] < start ?
                            delta_array[i] *= mOptions.mInc : delta_array[i] *= mOptions.mDec;
                    new_index = min_result == result_direction[i] ?
                            new_index = index[i] : new_index = new_index;
                }
                if (start > min_result) {
                    start = min_result;
                    int right_index = abs(new_index) - 1;
                    if (new_index < 0) {
                        starting_point[right_index] = std::max(starting_point[right_index] - delta_array[right_index],
                                                               left_border[right_index]);
                    } else {
                        starting_point[right_index] = std::min(starting_point[right_index] + delta_array[right_index],
                                                               right_border[right_index]);
                    }
                    delta_array[right_index] *= mOptions.mInc;
                }
            }
            return start;
        }
    };
}

#endif //PANTHER_MASTER3_PARALLCOORDDESCMETHOD_H
