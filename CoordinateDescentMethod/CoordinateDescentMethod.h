#include <sstream>
#include <vector>
#include <functional>
#include <memory>
#include <common/bbsolver.hpp>
#include <common/vec.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iomanip>
#pragma once

namespace panther {
    template <typename T> class CoordinateDescentMethod : public BlackBoxSolver<T> {
    public:
        struct Options {
            T delta = 1e-1;
            T mMinStep = 1e-3;
            T mDec = 0.5;
            T mInc = 1.5;
        } mOptions;

        bool check_border(int N, const T* left_border, const T* right_border) {
            //left >= right
            for (int i = 0; i < N; ++i) {
                if (left_border[i] > right_border[i]) {
                    return false;
                }
            }
            return true;
        }

        bool is_in_box(int n, const T* vec_for_check, const T* left_border, const T* right_border) {
            for (int i = 0; i < n; i++) {
                if (vec_for_check[i] > right_border[i] || vec_for_check[i] < left_border[i]) {
                    return false;
                }
            }
            return true;
        }

        //покоординатный спуск
        T search(int N, T* start, const T* left_border, const T* right_border, const std::function<T(const T * const)> &f) override {
            //проверка коррректности введенных данных:
            if (!check_border(N, left_border, right_border)) {
                throw std::invalid_argument("Границы интервала введены некорректно\n");
            }
            if (!is_in_box(N, start, left_border, right_border)) {
                throw std::invalid_argument("Точка не принадлежит интервалу\n");
            }
            //данные корректны

            std::vector<T> deltas_array(N, mOptions.delta);
            auto maxStep = [&deltas_array, N]() {
                T max_delta = 0;
                for (int i = 0; i < N; i++) {
                    max_delta = std::max(max_delta, deltas_array[i]);
                }
                return max_delta;
            };

            T present_value, previous_value, save_start;
            present_value = f(start);

            while (maxStep() > mOptions.mMinStep) {// уменьшаем погрешность, пока это возможно
                for (int i = 0; i < N; ++i) {    // идём по всем направлениям
                    previous_value = start[i]; //сохраняем изначальное значение start[i]
                    start[i] += deltas_array[i]; //идем вдоль оси
                    if (start[i] > right_border[i]) { //если мы сделали шаг вправо и вышли за правую границу
                        start[i] = right_border[i];
                    }
                    T right_d = f(start); //по дефолту скажем, что правильное направление вдоль оси
                    if (right_d >= present_value) { //если вдоль оси не уменьшается значение
                        start[i] -= 2 * deltas_array[i]; //пробуем пойти против оси
                        if (start[i] < left_border[i]) { //если сделали шаг влево и вышли за левую границу
                            start[i] = left_border[i];
                        }
                        T against_d = f(start); //значение функции против оси
                        if (against_d >= present_value) { //если против оси не уменьшается значение
                            //уменьшаем дельту, если никуда не пошли
                            deltas_array[i] *= mOptions.mDec;
                            start[i] = previous_value; //возвращаем изначальный старт
                            continue;
                        }
                        right_d = against_d;
                        deltas_array[i] *= -1; //умножаем дельту на -1, тк идем против оси
                    }
                    save_start = start[i];
                    while (right_d < present_value) { //пока значение улучшается, идем в выбранном направлении
                        deltas_array[i] *= mOptions.mInc;
                        present_value = right_d;
                        save_start = start[i];
                        start[i] += deltas_array[i];
                        if (start[i] < left_border[i] || start[i] > right_border[i]) { //если вышли за границы, то выходим из цикла
                            break;
                        }
                        right_d = f(start);
                    }
                    start[i] = save_start;
                    deltas_array[i] = deltas_array[i] < 0 ? deltas_array[i] * -1 : deltas_array[i]; //возвращаем знак дельты
                }
            }
            return present_value;
        }
    };
}
