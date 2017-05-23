#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <string.h>
#include "mersenne.cpp"

// #define CHECK_WRITE

#define energy(b) (2*(b)-3.0*(Lx*Ly*Lz))

int Lx, Ly, Lz;

double Neighbours( double ***s, int i, int j, int k )
{
    double res;

    if(i==0) res = s[Lx-1][j][k];
    else res = s[i-1][j][k];

    if(i==Lx-1) res = res + s[0][j][k];
    else res = res + s[i+1][j][k];

    if(j==0) res = res + s[i][Ly-1][k];
    else res = res + s[i][j-1][k];

    if(j==Ly-1) res = res + s[i][0][k];
    else res = res + s[i][j+1][k];

    if(k==0) res = res + s[i][j][Lz-1];
    else res = res + s[i][j][k-1];

    if(k==Lz-1) res = res + s[i][j][0];
    else res = res + s[i][j][k+1];

    return res;
}

int main()
{
    int i, j, k, N, a, b=0, b_new,
    n, skip, min_steps, mcs, rng, top_b, *hist, it_count;

    double *m, *m2, *m_sum, *m2_sum, m_sum_no_arr, m2_sum_no_arr;

    double f, fmin, dec_pow, flat_thres,
    r, ksi1, ksi2, ksi3, p,
    *g, ***s_x, ***s_y, ***s_z,
    Delta[36]={0.01,
        0.522701295579425,
        0.581317791984826,
        0.63601096704108,
        0.686780820748186,
        0.733627353106146,
        0.776550564114959,
        0.815550453774625,
        0.850627022085143,
        0.881780269046515,
        0.90901019465874,
        0.932316798921817,
        0.951700081835748,
        0.967160043400532,
        0.978696683616168,
        0.986310002482658,
        0.99,
        0.989766676168195,
        0.985610030987244,
        0.977530064457145,
        0.9655267765779,
        0.949600167349507,
        0.929750236771967,
        0.90597698484528,
        0.878280411569447,
        0.846660516944466,
        0.811117300970338,
        0.771650763647063,
        0.728260904974641,
        0.680947724953072,
        0.629711223582357,
        0.574551400862494,
        0.515468256793484,
        0.452461791375326,
        0.385532004608023,
        0.314678896491571},
    arcsin_Theta=1.0;

    it_count = 0;

    FILE *OUT;

    OUT = fopen("config.cfg", "r");
    fscanf(OUT, "%d", &Lx);
    fscanf(OUT, "%d", &Ly);
    fscanf(OUT, "%d", &Lz);
    fscanf(OUT, "%d", &rng);
    fscanf(OUT, "%lf", &f);
    fscanf(OUT, "%lf", &fmin);
    fscanf(OUT, "%lf", &dec_pow);
    fscanf(OUT, "%lf", &flat_thres);
    fscanf(OUT, "%d", &min_steps);
    fscanf(OUT, "%d", &skip);
    fclose(OUT);

    N = Lx*Ly*Lz;

    /** UNCOMMENT FOR ANISOTROPY **/
    //arcsin_Theta = asin(1.0-Delta[Lz]);

    if(!(Lz%2))
        top_b=Lx*Ly*Lz+1;
    else
        top_b=Lx*Ly*(Lz-1)+1;//if that happens then everything is bad

    hist = new int [top_b];
    g = new double [top_b];

    m = new double [top_b];
    m_sum = new double [top_b];

    m2 = new double [top_b];
    m2_sum = new double [top_b];

    s_x = new double** [Lx];
    s_y = new double** [Lx];
    s_z = new double** [Lx];

    for(i=0;i<Lx;i++)
        {
            s_x[i] = new double* [Ly];
            s_y[i] = new double* [Ly];
            s_z[i] = new double* [Ly];

                for(j=0;j<Ly;j++)
                    {
                        s_x[i][j] = new double [Lz];
                        s_y[i][j] = new double [Lz];
                        s_z[i][j] = new double [Lz];
                    }
        }

    for(i=0;i<Lx;i++)
        {
            for(j=0;j<Ly;j++)
                {
                    for(k=0;k<Lz;k++)
                        {
                            s_x[i][j][k] = 0.0;
                            s_y[i][j][k] = 0.0;
                            s_z[i][j][k] = 1.0;
                        }
                }
        }

    for(i=0;i<top_b;i++)
        g[i]=1.0;

    CRandomMersenne Mersenne(time(0));

    #ifdef CHECK_WRITE
        goto m1;
    #endif

    while(f>fmin)
        {
            printf("f=%.9lf\n", f);

            double lnf=log(f);
            int c, cont=1;

            for(a = 0;a < top_b; a++)   {
                hist[a] = 0;
            }

            n=0;
            c=skip+1;
            mcs=0;

            do
                {
                    if(mcs%100000==0)
                        printf("mcs %d\n", mcs);

                    for(a=0;a<N;a++)
                        {
                            i = Mersenne.IRandomX(0, Lx-1);
                            j = Mersenne.IRandomX(0, Ly-1);
                            k = Mersenne.IRandomX(0, Lz-1);

                            b_new=b+(int)((s_x[i][j][k] * Neighbours(s_x,i,j,k)) +
                            (s_y[i][j][k] * Neighbours(s_y,i,j,k)) +
                            (s_z[i][j][k] * Neighbours(s_z,i,j,k)));

                           if(b_new >= 0 && b_new < top_b )
                                {
                                    p = exp(g[b] - g[b_new]);

                                    if((p>=1.0)||(Mersenne.Random() < p))
                                        {
                                            do
                                                {
                                                    ksi1 = (double) Mersenne.IRandomX(-rng,rng)/rng;
                                                    ksi2 = (double) Mersenne.IRandomX(-rng,rng)/rng;
                                                    ksi3 = (double) Mersenne.IRandomX(-rng,rng)/rng;

                                                    r = sqrt(ksi1*ksi1 + ksi2*ksi2 + ksi3*ksi3);

                                                } while(r>1.0);

                                            s_x[i][j][k] = ksi1*arcsin_Theta/r;
                                            s_y[i][j][k] = ksi2*arcsin_Theta/r;
                                            s_z[i][j][k] = ksi3/r;

                                            b = b_new;
                                        }

                                    g[b]+=lnf;

                                    m[b]+=s_x[i][j][k];
                                    m[b]+=s_y[i][j][k];
                                    m[b]+=s_z[i][j][k];

                                    m2[b]+=s_x[i][j][k]*s_x[i][j][k];
                                    m2[b]+=s_y[i][j][k]*s_y[i][j][k];
                                    m2[b]+=s_z[i][j][k]*s_z[i][j][k];

                                    hist[b] += 1;
                                    n++;
                                }
                        }

                    mcs++;
                    c++;

                    if((mcs>=min_steps) && (c>=skip))
                        {
                            int div;
                            c=0;
                            div=top_b>Lx*Ly*Lz-1 ? top_b-2 : top_b-1;

                            for(a=0;
                                (a<top_b) && ((a==1) || (a==Lx*Lx*Lz-1) || ((double)hist[a]/n*div > flat_thres));
                                    a++);
                            if(a==top_b)
                                cont=0;
                            else
                                cont=1;
                        }
                    else
                        cont=1;

                    it_count++;

                    if(mcs >= 300000)   {

                        break;

                    }

                } while(cont);

                /** alternative flatness check **/
                    /*if((mcs >= min_steps) && (c >= skip))
                        {
                            c = 0;

                            int h_count = 0;
                            double h_delt;
                            double h_sum = 0.0;

                            for(int i = 0; i < top_b; i++)
                                {
                                    if(hist[i] != 0)
                                        {
                                            h_sum += (double)hist[i]/min_steps;
                                            h_count++;
                                        }

                                }

                            cont = 0;

                            for (int i = 0; i < top_b; i++)
                                {
                                    if(hist[i] != 0)
                                        {
                                            h_delt = (double)hist[i]/(h_sum/h_count)/min_steps;
                                            if((h_delt <= flat_thres) || (h_delt >= 1.0+flat_thres))
                                                cont = 1;
                                        }
                                }
                        }
                    else
                        cont = 1;
                }while(cont);*/

            std::cout << "top_b = " << top_b << std::endl;

            double min_gi = g[0];

            for(int i = 0; i < top_b; i++)  {

                if(min_gi > g[i])   {

                    min_gi = g[i];

                }

            }

            for(int i = 1; i < top_b; i++)    {
    
                g[i] -= min_gi;
    
                m_sum[i] += m[i];
                // m_sum_no_arr += m[i];
    
                // std::cout << m_sum[i] << std::endl;
    
                m[i] = 0.0;
    
                m2_sum[i] += m2[i];
    
                // m2_sum_no_arr += m2[i];
    
                std::cout << m2_sum[i] << std::endl;
    
                m2[i] = 0.0;
    
                g[0] = 0.0;
    
            }

            f = pow(f,dec_pow);
        }

    OUT = fopen("data.dat", "w+");
        {
            m1:

            for(a=0;a<top_b;a++)
                {
                    if((a!=1)&&(a!=Lx*Ly*Lz-1))
                        {
                            fprintf(OUT, "%d", a-top_b/2);
                            fprintf(OUT, "\t%.9lf", g[a]+log(100000.0));
                            fprintf(OUT, "\t%d\n", hist[a]);
                        }
                }

                    char * filename_out_td = new char [100];  // termodinamical
                    char * filename_out_ds = new char [100];  // density of states
                    char * filename_out_mg = new char [100];  // magnet

                    double EE, EE2, GE, Ut, Ft, St, Ct, lambdatemp, lambda, MM, MM2, Mt, Xt;

                    std::ofstream out_f_td, out_f_ds, out_f_mg, time_f;

                    std::stringstream ss;

                    ss.str("");
                    ss << "results/TermodinamicalStat_HEIS_3D_L=" << Lx << ".dat";
                
                    strcpy(filename_out_td ,ss.str().c_str());
                
                    out_f_td.open(filename_out_td, std::ios::out);
                    if(!out_f_td) std::cout << "Cannot open " << filename_out_td << ".Check permissions or free space! \n";
                    out_f_td << "T\tUt\tFt\tSt\tCt\n";
                
                    ss.str("");
                    ss << "results/DensityStat_HEIS_3D_L=" << Lx << ".dat";
                
                    strcpy(filename_out_ds, ss.str().c_str());
                
                    out_f_ds.open(filename_out_ds, std::ios::out);
                    if(!out_f_ds) std::cout << "Cannot open " << filename_out_ds << ".Check permissions or free space! \n";
                    out_f_ds << "i\tE(i)\tg[i]\thist[i]\n";

                    ss.str("");
                    ss << "results/MagnetStat_HEIS_3D_L=" << Lx << ".dat";
                
                    strcpy(filename_out_mg ,ss.str().c_str());
                
                    out_f_mg.open(filename_out_mg, std::ios::out);
                    if(!out_f_mg) std::cout << "Cannot open " << filename_out_mg << ".Check permissions or free space! \n";
                    
                    std::cout << "LX LY LZ N top_b: " << Lx << " : " << Ly << " : " << Lz << " : " << N << " : " << top_b << "\n";

                    for(int i = 0; i < top_b; i++)  {

                        std::cout << "M[i] = " << m[i] << std::endl;

                    }

                    for(double T = 0; T <= 8; T += 0.01)  {
                
                        EE = 0;
                        EE2 = 0;
                        GE = 0;

                        MM = 0;
                        MM2 = 0;
                
                        lambda = 0;
                        lambdatemp = 0;
                
                        for(int i = 0; i < top_b; i++)  {
                            // if((i!=N-1) && (hist[i]!=0))    {

                                lambdatemp = g[i] - energy(i)/T;
                                if(lambdatemp > lambda) lambda = lambdatemp;

                            // }
                        }
                
                        for(int i = 0; i < top_b; i++) {
                            // if((i!=N-1) && (hist[i]!=0))    {
                                
                                EE  += energy(i)*exp(g[i]-(energy(i))/T-lambda);
                                EE2 += energy(i)*energy(i)*exp(g[i]-(energy(i))/T-lambda);
                                GE  += exp(g[i]-energy(i)/T-lambda);

                                MM  += m_sum[i]*exp(g[i]-(energy(i))/T-lambda);

                                MM2 += m2_sum[i]*exp(g[i]-(energy(i))/T-lambda);

                                // std::cout << "GE = " << GE << std::endl;

                            // }
                        }

                        std::cout << "MM = " << MM << std::endl;
                        std::cout << "MM2 = " << MM2 << std::endl;
                        
                        // MM2 = MM2/GE;
                
                        Ut = EE/GE;
                        Ft = -T*lambda-(T)*log(GE);
                        St = (Ut-Ft)/T;
                        Ct = ((EE2/GE)-Ut*Ut)/(T*T);

                        Mt = MM/GE;
                        Xt = (Mt*Mt-(MM2/GE))/T;
                
                        // if((Ut==Ut)==true&&(Ft==Ft)==true&&(St==St)==true&&(Ct==Ct)==true) // Не пишем NaN 
                        out_f_td << std::fixed << T << "\t" << Ut/N << "\t" << Ft/N << "\t" << St/N << "\t" << Ct/N << "\n";
                        out_f_mg << std::fixed << T << "\t" << Mt/N << "\t" << Xt/N << "\n";

                    }

                    for(int i = 0; i < top_b; i++)  {

                        out_f_ds << std::fixed << std::setprecision(6) << i << "\t" << energy(i) << "\t" << g[i] << "\t" << hist[i] << "\n";

                    }   

                    out_f_td.close();
                    out_f_ds.close();
                    out_f_mg.close();

        }
    fclose(OUT);

    OUT = fopen("temp.cfg","w+");
    fprintf(OUT, "%d\n", Lx);
    fprintf(OUT, "%d\n", Ly);
    fprintf(OUT, "%d\n", Lz);
    fprintf(OUT, "%d\n", top_b);
    fclose(OUT);

    printf("\a \a");
    return 0;
}
