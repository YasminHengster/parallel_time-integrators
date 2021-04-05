/* Author: Tianming Bai
    Based on Hossein's python code.
    Not well commented.
*/
#include <cmath>
#include <stdio.h>
#include "mpi.h"
#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>

const int N_posi = 200, N_nega = 200;
const double q_posi = 1/(double)N_posi, q_nega = -1/(double)N_nega;
const double m_posi = (double)1000/(double)N_posi, m_nega = 1/(double)N_nega;
const double d = 0.05;


double fx_posi(int i, double y[2*(N_posi+N_nega)]){
    //the RHS of the system of ODEs y'=f(t,y)
    double ff = y[2*i+1];
    return ff;
}

double fv_posi(int ii, double y[2*(N_posi+N_nega)]){
    //the RHS of the system of ODEs y'=f(t,y)
    double ff = 0, temp1 = 0, temp2 = 0;
    for(int k = 0; k<N_posi;k++){
        temp1 = temp1 + (y[2*ii]-y[2*k])/sqrt((y[2*ii]-y[2*k])*(y[2*ii]-y[2*k])+d*d);

    }
    for(int k = 0; k<N_nega;k++){
        temp2 = temp2 + (y[2*ii]-y[2*N_posi+2*k])/sqrt((y[2*ii]-y[2*N_posi+2*k])*(y[2*ii]-y[2*N_posi+2*k])+d*d);
    }
    
    ff = q_posi * (q_posi*temp1+q_nega*temp2)/ m_posi;
    
    return ff;
}

double fx_nega(int i, double y[2*(N_posi+N_nega)]){
    //the RHS of the system of ODEs y'=f(t,y)
    double ff = y[2*N_posi+2*i+1];
    return ff;
}

double fv_nega(int i, double y[2*(N_posi+N_nega)]){
    //the RHS of the system of ODEs y'=f(t,y)
    double ff = 0, temp1 = 0, temp2 = 0;
    for(int k = 0; k<N_posi;k++){
        temp1 = temp1 + (y[2*N_posi+2*i]-y[2*k])/sqrt((y[2*N_posi+2*i]-y[2*k])*(y[2*N_posi+2*i]-y[2*k])+d*d);
    }
    for(int k = 0; k<N_nega;k++){
        temp2 = temp2 + (y[2*N_posi+2*i]-y[2*N_posi+2*k])/sqrt((y[2*N_posi+2*i]-y[2*N_posi+2*k])*(y[2*N_posi+2*i]-y[2*N_posi+2*k])+d*d);
    }
    ff = q_nega * (q_posi*temp1 + q_nega*temp2) / m_posi;

    return ff;
}


int main(){
    /*    Inputs:
    ff: the RHS of the system of ODEs y'=f(t,y)
    T:  integration interval[0,T]
    y0: initial condition
    N:  number of nodes
    M: the number of points in calculating quadraure integral
    (and also the number of steps used in Adam-Bashforth predictor)
    or number of correction loops PLUS the prection loop

    Output:
    t: time vector
    yy: solution as a function of time
    */

    
    double x_i;
    double T = 100, y0[2*(N_posi+N_nega)];
    int N = 100, M = 4;

    for(int i=0; i<N_posi;i++){
        y0[2*i] = (i-0.5)/N_posi;
        y0[2*i+1] = 0;
    }
    for(int i=0;i<N_nega;i++){
        x_i = (i-0.5)/N_nega;
        y0[2*N_posi+2*i] = x_i;
        y0[2*N_posi+2*i+1] = sin(6*M_PI *x_i);
    }

    double h = T/N; //time step

    //initialize MPI    
    int nproc, rank;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
    double tstart, tend, tend2;
    double runtime, runtime2; 
    MPI_Barrier(MPI_COMM_WORLD);

    tstart = MPI_Wtime();


    double F0[2*(N_posi+N_nega)];
    for(int i = 0; i < N_posi; i++){
        F0[2*i] = fx_posi(i,y0);
        F0[2*i+1] = fv_posi(i,y0);

    }

    for(int i = 0; i < N_nega; i++){
        F0[2*N_posi+2*i] = fx_nega(i,y0);
        F0[2*N_posi+2*i+1] = fv_nega(i,y0);
    }

            
    /* ================== INITIAL PART (1) ==================
    # for this part the predictor and correctors step up to M points in time
    # ** predictor ** uses Runge-Kutta 4 */
    double Y_update[2*(N_posi+N_nega)];

    for(int iTime = 0; iTime<N; iTime++){
        for(int j = 0; j<2*(N_posi+N_nega);j++){
            y0[j] = y0[j] + F0[j]*h;
        }
        
  
        for(int j = 0; j<N_posi;j++){
            F0[2*j] = fx_posi(j,y0);
            F0[2*j+1] = fv_posi(j,y0);
        }
        for(int j = 0; j<N_nega;j++){
            F0[2*N_posi+2*j] = fx_nega(j,y0);
            F0[2*N_posi+2*j+1] = fv_nega(j,y0);
        }   

    }
   

    tend2 = MPI_Wtime();
    runtime = tend2 - tstart; 

    if(rank==0){        
        cout<<"t: "<<runtime << endl;
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            cout<<y0[j]<<' ';
        }
        cout<<endl;        
    }


    MPI_Finalize();

    if (rank == 0){          
        fstream myfile;
        string name("EU_100");
        string txt(".txt");
        name = name + txt;    
        myfile.open(name,fstream::out);
        for (int i = 0; i < 2*(N_posi+N_nega); i++) {
            myfile << fixed << setprecision(12)<<y0[i] << ",";
            myfile << endl; 
        }    
        myfile.close();
    } 


    return 0;
}