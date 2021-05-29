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
#include <thread>
#include <utility>
#include <mutex>
#include <future>

#ifndef PANTHER_MASTER_OLDPARALLCOORDDESCMETHOD_H
#define PANTHER_MASTER_OLDPARALLCOORDDESCMETHOD_H

namespace panther {
    template <typename T> class OldparallCoordinateDescentMethod : public BlackBoxSolver<T> {
    public:
        struct Options {
            T delta = 1e-1;
            T mMinStep = 1e-2;
            T mDec = 0.5;
            T mInc = 1.5;
        } mOptions;

        bool check_border(int N, const T* left_border, const T* right_border) {
            for (int i = 0; i < N; ++i) {
                if (left_border[i] > right_border[i]) {
                    return false;
                }
            }
            return true;
        }

        bool is_in_box(int n, const T* vec_for_check, const T* left_border, const T* right_border) {
            for (int i = 0; i < n; i++) {
                if (vec_for_check[i] >= right_border[i] || vec_for_check[i] <= left_border[i]) {
                    return false;
                }
            }
            return true;
        }

        void search_one_direction(const T left_border, const T right_border, T & deltas_array,int i,
                const std::function<T(const T * const)> &f, T* start, std::vector<T>& result, std::vector<T>& answer,
                int N) {

            T present_value;
            T new_start[N];
            for (int j = 0; j < N; ++j) {
                new_start[j] = start[j];
            }
            present_value = f(new_start);

            new_start[i] = std::min(new_start[i] + deltas_array, right_border);

            T right_d = f(new_start); //по дефолту скажем, что правильное направление вдоль оси

            if (right_d >= present_value) { //если вдоль оси не уменьшается значение
                new_start[i] = std::max(new_start[i] - 2 * deltas_array, left_border);

                right_d = f(new_start); //значение функции против оси
                if (right_d >= present_value) { //если против оси не уменьшается значение
                    //уменьшаем дельту, если никуда не пошли
                    deltas_array *= mOptions.mDec;
                    return;
                }
                deltas_array *= -1; //умножаем дельту на -1, тк идем против оси
            }
            T save_start = new_start[i];
            while (right_d < present_value) { //пока значение улучшается, идем в выбранном направлении
                deltas_array *= mOptions.mInc;
                present_value = right_d;
                save_start = new_start[i];
                new_start[i] += deltas_array;
                if (new_start[i] + deltas_array < left_border ||
                new_start[i] + deltas_array > right_border) { //если вышли за границы, то выходим из цикла
                    break;
                }
                right_d = f(new_start);
            }
            result[i] = save_start;
            deltas_array = deltas_array < 0 ? deltas_array * -1 : deltas_array;//возвращаем знак дельты
            answer[i] = present_value;
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
                    max_delta = std::max((max_delta), (deltas_array[i]));
                }
                return max_delta;
            };
            std::vector <std::thread> threads;
            std::vector<T> result(N);
            std::vector<T> answer(N);
            for (int i = 0; i < N; ++i) {
                result[i] = start[i];
            }
            T help = f(start);
            while (maxStep() > mOptions.mMinStep) {// уменьшаем погрешность, пока это возможно
                for (int i = 0; i < N; ++i) {// идём по всем направлениям
                    std::thread t(std::thread(&OldparallCoordinateDescentMethod::search_one_direction, this, left_border[i],
                            right_border[i], std::ref(deltas_array[i]), i, f, start, std::ref(result), std::ref(answer),
                            N));
                    threads.push_back(std::move(t));
                }
                for(std::thread &t: threads){
                    if (t.joinable()) {
                        t.join();
                    }
                }
                std::pair<T, int> pair_ans = {answer[0], 0};
                for (int i = 0; i < N; ++i) {
                    if (answer[i] < pair_ans.first) {
                        pair_ans.first = answer[i];
                        pair_ans.second = i;
                    }
                }
                start[pair_ans.second] = result[pair_ans.second];
                help = answer[pair_ans.second];
            }
            return help;
        }
    };
}

#endif //PANTHER_MASTER_OLDPARALLCOORDDESCMETHOD_H