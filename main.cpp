#include <list>
#include <algorithm>
#include <dirent.h>
#include <omp.h>
#include "landau-wang.h"

using namespace std;

int main(int argc, char *argv[])  {

    long omp_time_start, omp_time_end, runtime_sec;
    stringstream command, time_filename;

    int L, MIN_STEPS, IT_STEPS, MAGNET_IT_AFTER;
    double F, F_MIN, F_THRESHOLD, T_MIN, T_MAX;
    string SMM;

    bool uncorrect_parameters = false;

    DIR *path_to_files_with_parameters;
    struct dirent *rec_of_file;

    int count_of_input_files, count_of_parameters;
    std::string current_parameter;
    std::ifstream list_of_parameters;

    list<string> list_of_input_files;
    list<string>::iterator input_files_iterator;

    count_of_input_files = 0;
    count_of_parameters  = 0;

    stringstream buffer;
    string filename;

    if((path_to_files_with_parameters = opendir("./input_parameters/")) != NULL)   {

        while((rec_of_file = readdir(path_to_files_with_parameters)) != NULL)   {

            count_of_input_files++;

            buffer.str("");
            buffer << rec_of_file->d_name;

            list_of_input_files.insert(list_of_input_files.begin(), buffer.str().c_str());

        }

    }

    bool increment_iterator = true;
    // Пройдем весь двусвязный список от начала до конца с помощью итератора
    for(auto it = list_of_input_files.begin(); it != list_of_input_files.end(); it++)  {
        // Делаем инкремент итератора +2 в начале для того чтобы исключить "." и ".." из списка имен файлов
        if(increment_iterator) {

            increment_iterator = false;
            // Нужно для перемещения итератора на определенную позицию списка
            advance(it,0);

        }

        count_of_input_files++;

        buffer.str("");
        buffer << "input_parameters/" << (string)*it;

        filename = buffer.str();
        replace(filename.begin(), filename.end(), '\n', ' ');

        // cout << filename << "\n";

        if(!(list_of_parameters.bad())) {

            list_of_parameters.close();

        }    

        list_of_parameters.open(filename);

        if(!(list_of_parameters.bad()))   {
            
            std::cout << "-------------------------------------- \n";
            std::cout << " ... READ FILE " << filename << "\n\n";
            // std::cout << "-------------------------------------- \n";

            if(filename == ".") 

            count_of_parameters = 0;
    
            while(getline(list_of_parameters, current_parameter))   {
    
                count_of_parameters++;

                if(count_of_parameters == 1) {
                    
                    // replace(current_parameter.begin(), current_parameter.end(), '\n', ' ');
                    // current_parameter.erase(current_parameter.size() - 1);
                    SMM = current_parameter;
                    // SMM - simulation model and method
                    std::cout << count_of_parameters << " :: SMM = " << SMM << "\n";                
    
                }
    
                if(count_of_parameters == 2) { 
    
                    L   = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: L = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 3) {
    
                    F = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: F = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 4) {
    
                    F_MIN = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: F_MIN = " << current_parameter << "\n";
    
                }
    
                if(count_of_parameters == 5) {
    
                    MIN_STEPS   = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: MIN_STEPS = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 6) {
    
                    IT_STEPS   = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: IT_STEPS = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 7) {
    
                    F_THRESHOLD = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: F_THRESHOLD = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 8) {
    
                    MAGNET_IT_AFTER = std::stoi(current_parameter);
                    // Итерации по намагниченности после количества общих итераций, указанных данным параметром
                    std::cout << count_of_parameters << " :: MAGNET_IT_AFTER = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 9) {
    
                    T_MIN = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_MIN = " << current_parameter << "\n";                
    
                }

                if(count_of_parameters == 10) {
    
                    T_MAX = std::stoi(current_parameter);
                    std::cout << count_of_parameters << " :: T_MAX = " << current_parameter << "\n";                
    
                }
    
            }
            
            cout << "\ncount_of_parameters = " << count_of_parameters << "\n\n";
            // cout << "-------------------------------------- \n";
            
            if(count_of_parameters < 7)    {

                cout << "One of the file have errors: \"" << filename << "\". Try to set correct count of parameters (10!). ABORT.\n\n";
                cout << "EXAMPLE:\n";
                cout << "-------------- \n";
                cout << "LW3D_HEIS\n8\n2.7182818284\n1.0000001\n10000\n0.8\n0\n1.5\n8\n";
                cout << "-------------- \n";
                cout << "DESCRIPTION:\n";
                cout << "-------------- \n";
                cout << "1 :: SMM                 (Simulation Method + Model e.g. LW3D_HEIS - [Landau-Wang + 3D Heisenberg model)\n";
                cout << "2 :: L                   (Size of the lattice)\n";
                cout << "3 :: F                   (Begin factor)\n";
                cout << "4 :: F_MIN               (End factor)\n";
                cout << "5 :: MIN_STEPS           (Minimal count of steps after which we check flatness\n";
                cout << "6 :: IT_STEPS            (Every IT_STEPS count of steps we check flatness\n";
                cout << "7 :: F_THRESHOLD         (Flat criteria)\n";
                cout << "8 :: MAGNET_IT_AFTER     (Count of steps after which we can estimate magnetisation)\n";
                cout << "9 :: T_MIN               (Begin of interval where steps for TEMP equal 0.1)\n";
                cout << "10 :: T_MAX              (End of interval where steps for TEMP equal 0.1)\n";
                std::cout << "-------------- \n";
                exit(1);

            }

            if(count_of_parameters == 10) {

                uncorrect_parameters = false;

                if((SMM != "LW3D_HEIS") && (SMM != "LW2D_ISING"))   {

                    uncorrect_parameters = true;
                    cout << "1 :: SMM = " << SMM << ". Set correct simulation method parameter! Must be [LW3D_HEIS, ... ]!\n";

                }

                if((L%2)!=0)    {
                    
                    uncorrect_parameters = true;
                    cout << "2 :: L = " << L << " .Set correct L! Must be devided on two!\n";
                    
                }

                if(F_MIN < 1.0)    {

                    uncorrect_parameters = true;
                    // В работе Ландау-Ванга было указано значение "f" равное экспоненте, а "f_min" около единицы, но больше
                    cout << "4,3 :: F_MIN = " << F_MIN << ", F_MIN usually near 1, but must be bigger! F = " << F << ", F usually equal exp()";

                }

                if(T_MIN < 0 || T_MAX < 0)   {

                    uncorrect_parameters = true;
                    cout << "8,9 :: T_min and T_max cannot be smaller than zero!";

                }

                if((T_MIN > T_MAX)) {

                    uncorrect_parameters = true;
                    cout << "8,9 :: T_min cannot be bigger than T_max!";

                }

                if(uncorrect_parameters) exit(1);
                else cout << " ... CONFIG LOOKING FINE! ... \n";
                
            }

            // Izing3D Model("Izing3D", L, MCS, STAT, PP);

            std::cout << SMM << ",\t" << L << ",\t" << F << ",\t" << F_MIN << ",\t" << MIN_STEPS << ",\t" << IT_STEPS << ",\t" << F_THRESHOLD << ",\t" << MAGNET_IT_AFTER << ",\t" <<  T_MIN << ",\t" <<  T_MAX << std::endl;

            LANDAU_WANG LW(SMM, L, F, F_MIN, MIN_STEPS, IT_STEPS, F_THRESHOLD, MAGNET_IT_AFTER, T_MIN, T_MAX);

            LW.LANDAU_WANG_2D_ISING();

            // Model.setMethod(SMM);

            // Model.setTemperatureRange(T_min, T_max, T_prec_min, T_prec_max);

            // Model.Start();

        }

    }

    return 0;

}