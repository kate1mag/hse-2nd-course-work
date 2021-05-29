#include <gtest/gtest.h>
#include "CoordinateDescentMethod/CoordinateDescentMethod.h"
#include "brute/bruteforce.hpp"
#include <vector>
#include "advcoordesc/advancedcoordescent.hpp"
#include "parallCoordDescMethod.h"
#include "CoordinateDescentMethod_old.h"
#include <cmath>
#include <ctime>
#include <iterator>
#include <chrono>

const int N = 10;
double x1[N];   //начальная точка поиска минимума
double x2[N];   //начальная точка поиска минимума
double x3[N];   //начальная точка поиска минимума
double k[N];    //вектор коэфф
double left[N]; //левая граница поиска
double right[N];//правая граница
double b = 0;   //смещение

void print(double* x, double f) {
    std::cout << "\nFOUND\t" << std::fixed << std::setprecision(7) << f << "\n";
    std::cout << "At\t";
    std::copy(x, x + N , std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n\n";
}

void fill() {
    std::fill(left, left + N, 1);
    std::fill(right, right + N, 2);
    std::fill(x1, x1 + N, 1.5);
    std::fill(x2, x2 + N, 1.5);
    std::fill(x3, x3 + N, 1.5);
    std::fill(k, k + N, 2);
}

class Function { //родительский класс функций
public:
    int N;
    std::vector<double> a;
    double b;
    Function(int n = 0, double* x = {}, double b1 = 0) : N(n) {
        for (int i = 0; i < n; ++i) {
            a.push_back(x[i]);
        }
        b = b1;
    }
};

class Line: public Function { // класс линейных функций
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += x[i] * a[i];
        }
        sleep(1);
        return v + b;
    }
};

class Parabola: public Function { // класс парабол
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += x[i] * x[i] * a[i] * pow(-1, i); // седло
        }
        sleep(1);
        return v + b;
    }
};

class Trig: public Function { // класс тригонометрических функций
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += cos(x[i]) * sin(x[i]) * a[i] ;
        }
        sleep(1);
        return v + b;
    }
};
class TestLine : public ::testing::Test {
protected:
    void SetUp()
    {
        fill();
        line = new Line(N, k, b);

        time_t start, end;
        time(&start);
        res = adv.search(N, x1, left, right, std::ref(*line));
        time(&end);
        double seconds = difftime(end, start);
        std::cout << "время последовательного\n";
        printf("The time: %f seconds\n", seconds);


        time_t start1, end1;
        time(&start1);
        res_adv_coord_d = adv_coord_d.search(N, x2, left, right, std::ref(*line)); // минимум рабочего метода
        time(&end1);
        double seconds1 = difftime(end1, start1);
        std::cout << "время параллельного1\n";
        printf("The time: %f seconds\n", seconds1);


        time_t start2, end2;
        time(&start2);
        res_adv_coord_d0 = adv_coord_d0.search(N, x3, left, right, std::ref(*line));
        time(&end2);
        double seconds2 = difftime(end2, start2);
        std::cout << "время параллельного2\n";
        printf("The time: %f seconds\n\n", seconds2);

    }
    void TearDown()
    {
        delete line;
    }
    panther::CoordinateDescentMethod<double> adv;
    panther::parallCoordinateDescentMethod<double> adv_coord_d;
    panther::OldparallCoordinateDescentMethod<double> adv_coord_d0;
    double res;
    double res_adv_coord_d;
    double res_adv_coord_d0;
    Line *line;
};

class TestParabola : public ::testing::Test {
protected:
    void SetUp()
    {
        fill();
        parabola = new Parabola(N, k, b);
        time_t start, end;
        time(&start);
        res1 = adv1.search(N, x1, left, right, std::ref(*parabola));
        time(&end);
        double seconds = difftime(end, start);
        std::cout << "время последовательного\n";
        printf("The time: %f seconds\n", seconds);

        time_t start1, end1;
        time(&start1);

        res_adv_coord_d1 = adv_coord_d1.search(N, x2, left, right, std::ref(*parabola));
        time(&end1);
        double seconds1 = difftime(end1, start1);
        std::cout << "время параллельного1\n";
        printf("The time: %f seconds\n", seconds1);


        time_t start2, end2;
        time(&start2);
        res_adv_coord_d10 = adv_coord_d10.search(N, x3, left, right, std::ref(*parabola));
        time(&end2);
        double seconds2 = difftime(end2, start2);
        std::cout << "время параллельного2\n";
        printf("The time: %f seconds\n\n", seconds2);
    }
    void TearDown()
    {
        delete parabola;
    }
    panther::CoordinateDescentMethod<double> adv1;
    panther::parallCoordinateDescentMethod<double> adv_coord_d1;
    panther::OldparallCoordinateDescentMethod<double> adv_coord_d10;
    Parabola *parabola;
    double res1;
    double res_adv_coord_d1;
    double res_adv_coord_d10;
};


class TestTrig : public ::testing::Test {
protected:
    void SetUp()
    {
        fill();
        trig = new Trig(N, k, b);
        time_t start, end;
        time(&start);
        res2 = adv2.search(N, x1, left, right, std::ref(*trig));
        time(&end);
        double seconds = difftime(end, start);
        std::cout << "время последовательного\n";
        printf("The time: %f seconds\n", seconds);

        time_t start1, end1;
        time(&start1);
        res_adv_coord_d2 = adv_coord_d2.search(N, x2, left, right, std::ref(*trig));
        time(&end1);
        double seconds1 = difftime(end1, start1);
        std::cout << "время параллельного1\n";
        printf("The time: %f seconds\n", seconds1);


        time_t start2, end2;
        time(&start2);
        res_adv_coord_d20 = adv_coord_d20.search(N, x3, left, right, std::ref(*trig));
        time(&end2);
        double seconds2 = difftime(end2, start2);
        std::cout << "время параллельного2\n";
        printf("The time: %f seconds\n\n", seconds2);
    }
    void TearDown()
    {
        delete trig;
    }
    panther::CoordinateDescentMethod<double> adv2;
    panther::parallCoordinateDescentMethod<double> adv_coord_d2;
    panther::OldparallCoordinateDescentMethod<double> adv_coord_d20;
    Trig *trig;
    double res2;
    double res_adv_coord_d2;
    double res_adv_coord_d20;
};

TEST_F(TestLine, test1) {
double result = abs(res_adv_coord_d - res);
std::cout << "последовательный\n";
print(x1, res);
std::cout << "параллельный1\n";
print(x2, res_adv_coord_d);
std::cout << "параллельный2\n";
print(x3, res_adv_coord_d0);
ASSERT_LE(result, 1e-1); // проверка с заданной точностью
}

TEST_F(TestParabola, test2) {
double result = abs(res_adv_coord_d1 - res1);
std::cout << "последовательный\n";
print(x1, res1);
std::cout << "параллельный1\n";
print(x2, res_adv_coord_d1);
std::cout << "параллельный2\n";
print(x3, res_adv_coord_d10);
ASSERT_LE(result, 1e-1);
}

TEST_F(TestTrig, test3) {
double result = abs(res_adv_coord_d2 - res2);
std::cout << "последовательный\n";
print(x1, res2);
std::cout << "параллельный1\n";
print(x2, res_adv_coord_d2);
std::cout << "параллельный2\n";
print(x3, res_adv_coord_d20);
ASSERT_LE(result, 1e-1);
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}