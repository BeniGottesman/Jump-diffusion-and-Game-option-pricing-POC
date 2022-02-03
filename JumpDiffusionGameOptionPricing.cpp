#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <deque>
#include <quadmath.h>

using namespace std;

// first form
typedef vector<double> v1d;
typedef vector<v1d> v2d;
//typedef vector<v2d> v3d;

static double K = 100.0;
static int n = 10; //tree depth
static double sigma = 0.7;
static double T = 1.0;
static double r = 0.0;
static double penalty = 10;
static double S0 = 105;
static double lambda = 0.01;//0.05;

double f(double x);
double g(double x);
double gameoptionValue(long double gamma, bool OnlyBS);

int main()
{
    double Solution = 0.0;

    //std::cout << "-10000 = "<< exp(-10000.0)<< " -1000 = "<< exp(-1000.0)<< " -100 = "<< exp(-100.0)<<"\n";
    //std::cout << "9= "<< exp(9.0)<< " 10= "<< exp(10.0)<< " 100 = "<< exp(100.0)<<"\n";
    Solution = gameoptionValue(0.0, true);
    std::cout << "B&S/Tree  Solution = "<< Solution <<"\n";

    long double gamma = -100.0;

    Solution = gameoptionValue(gamma, false);
    std::cout << "gamma = "<< gamma << " Solution = "<< Solution <<"\n";

    gamma = -1000.0;

    Solution = gameoptionValue(gamma, false);
    std::cout << "gamma = "<< gamma << " Solution = "<< Solution <<"\n";

    gamma = -10.0;

    Solution = gameoptionValue(gamma, false);
    std::cout << "gamma = "<< gamma << " Solution = "<< Solution <<"\n";

    for (int i=2; i<10; i=i+1)
    {
        gamma = (long double) i;
        Solution = gameoptionValue(gamma, false);
        std::cout << "gamma = "<< gamma << " Solution = "<< Solution <<"\n";
    }

    gamma = 50.0;
    Solution = gameoptionValue(gamma, false);
    std::cout << "gamma = "<< gamma << " Solution = "<< Solution <<"\n";
    return 0;
}

double gameoptionValue(long double gammap, bool OnlyBS)
{

    v2d * V;
    v2d * V_plus = new v2d (n+1, v1d(n+1));


    //Boundary condition
    //std::cout << "Boundary condition\n";
    double discount = exp(-r*T);
    for(int j = -n, k = n; j <= n; j+=2, k--)//j = xi
    {
        for(int m=n; m>=0; m--)//node m=n ==> n jumps
        {
            //double tmp = f(S0*exp(sigma*sqrt((T)/(double)n)*(double)j)*exp(log(1.0 + gamma)*(double)m - lambda*gamma*T)/discount);
            double tmp = f(S0*exp(sigma*sqrt((T)/(double)n)*(double)j) * exp(gammap*(double)m - lambda*(exp(gammap)-1.0)*T)/discount);
            (*V_plus) [k][m] = discount*tmp;
        }
    }


    //Backward Induction
    //std::cout << "Backward Induction\n";
    double P_jump = 0.0;
    //double P_jump2 = 0.0;

    if (OnlyBS == false)
    {
        P_jump = (1.0 - exp(0.0-lambda*(exp(gammap)-1.0)*T/n))
                /
                (exp(gammap-lambda*(exp(gammap)-1.0)*T/n) - exp(0.0-lambda*(exp(gammap)-1.0)*T/n));
        /*double u = 1.0/(gammap*gammap*lambda*T/n + 1.0);
        double d = 1.0/u;
        P_jump = (1.0 - d)/(u - d);*/
        /*double gamma = exp(gammap)-1.0;
        P_jump2 = (1.0 - exp(0.0-lambda*(gamma)*T/n))
                /
                (exp(log(1.0+gamma)-lambda*(gamma)*T/n) - exp(0.0-lambda*(gamma)*T/n));*/
    }

    //double gamma = exp(gammap)-1.0;
    //double mu = 2.0*lambda*gamma/sigma*sqrt(T/n);
    double P_BM_up = ( 1.0 - exp(-sigma*sqrt(T/n)))/(exp(sigma*sqrt(T/n)) - exp(-sigma*sqrt(T/n)));//UP
    //double P_BM_up = ( ((double)n/((double)n + mu*lambda*T)) - exp(-sigma*sqrt(T/n)))/(exp(sigma*sqrt(T/n)) - exp(-sigma*sqrt(T/n)));//UP

    for(int j = n-1; j>=0; j--)//depth
    {
        V = new v2d (j+1, v1d(j+1));
        discount = exp(-r*(double)j*(T/n));
        for(int i=-j, k = j; i<=j; i+=2, k--)//BS node i=sum(xi)
        {
            for(int m=j; m>=0; m--)//jump node m=j ==> j jumps
            {
                //double tmp_S_j = S0*exp(sigma*sqrt(T/n)*(double)i)*exp(log(1.0+gammap)*(double)m)/discount;
                double tmp_S_j = S0*exp(sigma*sqrt(T/n)*(double)i)*exp(gammap*(double)m - lambda*(exp(gammap)-1.0)*T/n*(double)j)/discount;
                double Y = discount*f(tmp_S_j);
                double X = discount*g(tmp_S_j);

                double EJ = 0.0;

                if(OnlyBS==true)
                {
                    //EJ = P_BM_up * V[j+1][k+1][0] + (1.0-P_BM_up) * V[j+1][k][0];
                    EJ = P_BM_up * (*V_plus)[k][0] + (1.0-P_BM_up) * (*V_plus)[k+1][0];
                    (*V)[k][0] = EJ; //EU + discount!!!!!!!
                    //(*V)[k][0] = max(Y, EJ); //American
                    //(*V)[k][0] = min (X, max(Y, EJ));
                    //V_plus[k+1][0].pop_back();
                }
                else
                {
                    //double v_tmp1, v_tmp2;
                    //v_tmp1 = (*V_plus)[k+1][m+1].pop_back();
                    //v_tmp2 = (double)(*V_plus)[k+1]][m+1].pop_back();

                    EJ = P_BM_up * ((*V_plus)[k][m+1]*P_jump + (*V_plus)[k][m]*(1.0-P_jump)) +
                            (1.0-P_BM_up) * ((*V_plus)[k+1][m+1]*P_jump + (*V_plus)[k+1][m]*(1.0-P_jump));
                    //EJ = P_BM_up * (V[j+1][k+1][m+1]*P_jump + V[j+1][k+1][m]*(1.0-P_jump)) +
                    //       (1.0-P_BM_up) * (V[j+1][k][m+1]*P_jump + V[j+1][k][m]*(1.0-P_jump));
                    (*V)[k][m] = EJ; //EU
                    //(*V)[k][m] = max (Y, EJ); //American
                    //(*V)[k][m] = min (X, max(Y, EJ));
                    (*V_plus)[k+1].pop_back();
                    //(*V_plus)[k].pop_back();
                }
            }//m
            //(*V_plus).pop_back();
        }//k
        V_plus->clear();
        //delete V_plus;
        V_plus = V;
    }//j

    double Solution = (*V)[0][0];

    delete V_plus;
    //delete V;

    return Solution;
}


double f(double x)
{
    return max(x-K,0.0);
}


double g(double x)
{
    return f(x) + penalty;
}
