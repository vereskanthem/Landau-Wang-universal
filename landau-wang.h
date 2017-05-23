/*
 * landau-wang.h
 *
 *  Created on: 24 мая 2014 г.
 *      Author: nlare
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>

#ifndef LANDAU_WANG_H_
#define LANDAU_WANG_H_

class LANDAU_WANG {

public:
    LANDAU_WANG(std::string, int, double, double, int, int, double, int, double, double);
    ~LANDAU_WANG();

    int    LANDAU_WANG_2D_ISING();
    int    LANDAU_WANG_3D_HEIS();
    int    NEIGHBOUR_ISING(int,int);
    double NEIGHBOUR_HEIS_3D(double***,int,int,int);

private:

    // SHARE VARS

    std::string SM;          // Simulation model

    int skip;                // Количество шагов с неизменной энергией
    int mcs;
    int min_steps;

    int it_magnet_count;    // После какой итерации считаем намагниченность
    int it_count;           // Подсчет количества итераций

    double f, f_min, ln_f;  /* "f" - начальный множитель для энергетических уровней, 
                             * "f_min" - минимальное значение множителя
                             * на каждый принятый шаг f_m = (f_m-1)^1/2
                             */

    int b, b_new, top_b;    // Текущий энергетический уровень и последующий, top_b - предельный энергетический уровень

    int L, N;               // Размер решетки

    int *hist, n;           // Массив гистограммы энергий, т.е количества посещений данного энергетического состояния
                            // n - суммарное количество посещений

    double *g;              /* Массив плотности состояний
                             * Изначально каждый элемент массива принимается равным единице
                             */

    double prob;             // Вероятность изменения энергетического уровня
    double flat_threshold;   // Порог "плоскости" гистограммы

    double T, T_min, T_max;

    double time_b, time_e;

    double EE, EE2, GE, Ut, Ft, St, Ct, lambdatemp, lambda;

    std::ofstream out_f_td, out_f_ds, plot_f, script_f, time_f;

    char * filename_out_td;
    char * filename_out_ds;

    std::stringstream buffer, ss;

    // ISING 2D

    int **spin;             // Массив спинов

    std::ofstream test_g_f, graph_g_f, graph_sh_g;

    // HEISENBERG 3D

    int Lx, Ly, Lz; // Размеры решетки в направлении x,y,z

    float arcsin_Theta, radius;

    double ***spin_x_heis, ***spin_y_heis, ***spin_z_heis, ksi1, ksi2, ksi3, dec_pow;

};



#endif /* LANDAU_WANG_H_ */
