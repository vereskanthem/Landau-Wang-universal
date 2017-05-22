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
    LANDAU_WANG(std::string, int, float, float, int, int, float, int, float, float);
    ~LANDAU_WANG();

    int LANDAU_WANG_2D_ISING();
    int LANDAU_WANG_3D_HEIS();
    int NEIGHBOUR_ISING(int,int);

private:

    int **spin;             // Массив спинов
    int L;                  // Размер решетки

    int *hist;              // Массив гистограммы энергий, т.е количества посещений данного энергетического состояния
    double *g;              /* Массив плотности состояний
                             * Изначально каждый элемент массива принимается равным единице
                             */
    int b, b_new, top_b;    // Текущий энергетический уровень и последующий, top_b - предельный энергетический уровень
    double f, f_min, ln_f;  /* "f" - начальный множитель для энергетических уровней, 
                             * "f_min" - минимальное значение множителя
                             * на каждый принятый шаг f_m = (f_m-1)^1/2
                             */
    int skip;                // Количество шагов с неизменной энергией
    int mcs;
    int min_steps;
    double prob;             // Вероятность изменения энергетического уровня
    double flat_threshold;   // Порог "плоскости" гистограммы

    double time_b, time_e;

    double T, T_min, T_max;

    int it_count;

    std::ofstream test_g_f, graph_g_f, graph_sh_g;

    double EE, EE2, GE, Ut, Ft, St, Ct, lambdatemp, lambda;

    std::ofstream out_f_td, out_f_ds, plot_f, script_f, time_f;

    char * filename_out_td;
    char * filename_out_ds;

};



#endif /* LANDAU_WANG_H_ */
