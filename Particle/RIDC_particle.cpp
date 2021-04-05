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

    int Mm = M-1; //the number of correctors

    //initialize MPI    
    int nproc, rank;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    /*  Generates beta coefficients for Adam-Bashforth integrating scheme
    These coefficients are stored in reversed compared to conventional
    Adam-Bashforth implementations (the first element of beta corresponds to
    earlier point in time).    */
    double beta_coeff[M]= {-9./24, 37./24, -59./24, 55./24};
    double beta_coeff2[M-1] = {5./12, -16/12, 23./12};
    /*if(M==3){
        float beta_coeff[3] = {5./12, -16/12, 23./12};
        float beta_coeff2[2] = {-1./2, 3/2};
    }else if(M==4){
        float beta_coeff[4] = {-9./24, 37./24, -59./24, 55./24};
        float beta_coeff2[3] = {5./12, -16/12, 23./12};
    }else if(M==5){
        float beta_coeff[5] = {251./720, -1274./720, 2616./720, -2774./720, 1901./720};
        float beta_coeff2[4] = {-9./24, 37./24, -59./24, 55./24};
    }else if(M==6){
        float beta_coeff[6] = {-475./720, 2877./720, -7298./720, 9982./720, -7923./720, 4277./720};
        float beta_coeff2[5] = {251./720, -1274./720, 2616./720, -2774./720, 1901./720};
    }*/

    /*  Generates quadrature weights, 
    whose elements are defined as integrals of Lagrange 
    interpolating polynomials.  */
    double S[Mm][Mm+1] = {{9./24, 19./24, -5./24, 1./24},
                             {-1./24, 13./24, 13./24, -1./24},
                             {1./24, -5./24, 19./24, 9./24}};
    /*if(M==3){       
        float S[Mm][Mm+1] = {{5./12, 8./12, -1/12},
                             {-1/12, 8./12, 5./12}};
    }else if(M==4){
        cout<<'b'<<endl;
        float S[Mm][Mm+1] = {{9./24, 19./24, -5./24, 1./24},
                             {-1./24, 13./24, 13./24, -1./24},
                             {1./24, -5./24, 19./24, 9./24}};
        cout<<S[0][0]<<endl;
    }else if(M==5){
        float S[Mm][Mm+1] = {{251./720, 646./720, 264./720, 106./720, -19./720},
                             {-19./720, 346./720, 456./720, 74./720, 11./720},
                             {11./720, 74./720, 456./720, 346./720, -19./720},
                             {-19./720, 106./720, 264./720, 646./720, 251./720}};
    }*/

    double Svec[Mm+1];
    for(int i = 0; i<Mm+1; i++){
        Svec[i] = S[Mm-1][i];
    }

    // the final answer will be stored in yy
    double yy[2*(N_posi+N_nega)];
    for(int i = 0; i < 2*(N_posi+N_nega); i++){
        yy[i] = y0[i];
    }

     
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



    // the time vector
    double t[N+1];
    for(int i = 0; i<N+1; i++){
        t[i] = i*h;
    }
    // extended time vector
    double t_ext[N+M+1];
    for(int i = 0; i<N+M+1; i++){
        t_ext[i] = i*h;
    }

    /* F vector and matrice:
    # the RHS of ODE is evaluated and stored in this vector and matrix:
    # F1 [M x M]: first index is the order (0=prection, 1=first correction)
    # second index is the time (iTime)
    # Note F1 could have been [M-1 x M] as the first two rows are equal to each
    # other BUT we designed it as a place holder for future parallelisation */
    double F1[2*(N_posi+N_nega)][Mm][M];
    for(int i = 0; i<Mm+1; i++){
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            F1[j][i][0] = F0[j];
        }        
    }
    

    double F2[2*(N_posi+N_nega)];
    for(int i = 0; i<2*(N_posi+N_nega); i++){
        F2[i] = F0[i];
    }

    //Y2 [M] new point derived in each level
    double Y2[2*(N_posi+N_nega)][M];
    for(int i = 0; i<M; i++){
        for(int j = 0; j<2*(N_posi+N_nega);j++){
            Y2[j][i] = y0[j];
        }        
    }
         
            
    /* ================== INITIAL PART (1) ==================
    # for this part the predictor and correctors step up to M points in time
    # ** predictor ** uses Runge-Kutta 4 */
    double KK1[2*(N_posi+N_nega)], KK2[2*(N_posi+N_nega)], KK3[2*(N_posi+N_nega)], KK4[2*(N_posi+N_nega)];
    double Y_update[2*(N_posi+N_nega)]; // Store temporary vector.
    for(int iTime = 0; iTime<M-1; iTime++){
        // KK1
        for(int j = 0; j<2*(N_posi+N_nega);j++){
            KK1[j] = F1[j][0][iTime];
            Y_update[j] = Y2[j][0] + KK1[j]*h/2;
        }
        // KK2
        for(int j = 0; j<N_posi;j++){
            KK2[2*j] = fx_posi(j,Y_update);
            KK2[2*j+1] = fv_posi(j,Y_update);          
        }
        for(int j = 0; j<N_nega;j++){
            KK2[2*N_posi+2*j] = fx_nega(j,Y_update);
            KK2[2*N_posi+2*j+1] = fv_nega(j,Y_update);
        }

        for(int j = 0; j<N_posi;j++){
            Y_update[2*j] = Y2[2*j][0] + KK2[2*j]*h/2;
            Y_update[2*j+1] = Y2[2*j+1][0] + KK2[2*j+1]*h/2;  
        }
        for(int j = 0; j<N_nega;j++){
            Y_update[2*N_posi+2*j] = Y2[2*N_posi+2*j][0] + KK2[2*N_posi+2*j]*h/2;
            Y_update[2*N_posi+2*j+1] = Y2[2*N_posi+2*j+1][0] + KK2[2*N_posi+2*j+1]*h/2;            
        }

        // KK3
        for(int j = 0; j<N_posi;j++){
            KK3[2*j] = fx_posi(j,Y_update);
            KK3[2*j+1] = fv_posi(j,Y_update);
        }
        for(int j = 0; j<N_nega;j++){
            KK3[2*N_posi+2*j] = fx_nega(j,Y_update);
            KK3[2*N_posi+2*j+1] = fv_nega(j,Y_update);
        } 
        for(int j = 0; j<N_posi;j++){
            Y_update[2*j] = Y2[2*j][0] + KK3[2*j]*h;
            Y_update[2*j+1] = Y2[2*j+1][0] + KK3[2*j+1]*h;            
        }
        for(int j = 0; j<N_nega;j++){
            Y_update[2*N_posi+2*j] = Y2[2*N_posi+2*j][0] + KK3[2*N_posi+2*j]*h;
            Y_update[2*N_posi+2*j+1] = Y2[2*N_posi+2*j+1][0] + KK3[2*N_posi+2*j+1]*h;            
        }          
        //KK4
        for(int j = 0; j<N_posi;j++){
            KK4[2*j] = fx_posi(j,Y_update);
            KK4[2*j+1] = fv_posi(j,Y_update);
        }
        for(int j = 0; j<N_nega;j++){
            KK4[2*N_posi+2*j] = fx_nega(j,Y_update);
            KK4[2*N_posi+2*j+1] = fv_nega(j,Y_update);
        } 
        //Y2[][0]
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            Y2[j][0] = Y2[j][0] + h*(KK1[j] + 2*KK2[j] + 2*KK3[j] + KK4[j])/6;        
        }     
        //F1[][0][iTime]
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            Y_update[j] = Y2[j][0];        
        }  
        for(int j = 0; j<N_posi;j++){
            F1[2*j][0][iTime+1] = fx_posi(j,Y_update);
            F1[2*j+1][0][iTime+1] = fv_posi(j,Y_update);
        }
        for(int j = 0; j<N_nega;j++){
            F1[2*N_posi+2*j][0][iTime+1] = fx_nega(j,Y_update);
            F1[2*N_posi+2*j+1][0][iTime+1] = fv_nega(j,Y_update);
        }   

    }

    //** correctors ** use Integral Deffered Correction
    int ll; double temp;
    for(int iCor =  1; iCor<M-1; iCor++){
        ll = iCor - 1;
        for(int iTime = 0; iTime<M-1; iTime++){
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                temp = 0;
                for(int i = 0; i<M; i++){
                    temp = temp + S[iTime][i]*F1[j][ll][i];
                }  
                Y2[j][iCor] = Y2[j][iCor]+ h*(F1[j][iCor][iTime]-F1[j][ll][iTime]) + h * temp;  
            }
  

            for(int j = 0; j<2*(N_posi+N_nega); j++){
                Y_update[j] = Y2[j][iCor];        
            }  

            for(int j = 0; j<N_posi;j++){
                F1[2*j][iCor][iTime+1] = fx_posi(j,Y_update);
                F1[2*j+1][iCor][iTime+1] = fv_posi(j,Y_update);
            }
            for(int j = 0; j<N_nega;j++){
                F1[2*N_posi+2*j][iCor][iTime+1] = fx_nega(j,Y_update);
                F1[2*N_posi+2*j+1][iCor][iTime+1] = fv_nega(j,Y_update);
            }   
        }
    }

    // treat the last correction loop a little different
    for(int iTime = 0; iTime<M-1; iTime++){
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            temp = 0;
            for(int i = 0; i<M; i++){
                temp = temp + (S[iTime][i]*F1[j][M-2][i]);
            }
            Y2[j][M-1] = Y2[j][M-1] + h*(F2[j]-F1[j][M-2][iTime]) + h * temp;          
        }

        //Use Y2 update F2
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            Y_update[j] = Y2[j][M-1];        
        }  
        for(int j = 0; j<N_posi;j++){
            F2[2*j] = fx_posi(j,Y_update);
            F2[2*j+1] = fv_posi(j,Y_update);
        }
        for(int j = 0; j<N_nega;j++){
            F2[2*N_posi+2*j] = fx_nega(j,Y_update);
            F2[2*N_posi+2*j+1] = fv_nega(j,Y_update);
        } 
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            yy[j] = Y2[j][M-1];       
        }   
        
    }

    //================== INITIAL PART (2) ==================
    int iStep, iCor;
    for(int iTime = M-1; iTime<2*M-2; iTime++){
        iStep = iTime - (M-1);
        // prediction loop
        for(int j = 0; j<2*(N_posi+N_nega); j++){            
            temp = 0;
            for(int i = 0; i<M; i++){
                temp = temp + beta_coeff[i]*F1[j][0][i];
            }
            Y2[j][0] = Y2[j][0] + h*temp;            
        }
        //correction loop
        for(int ll = 0; ll<iStep; ll++){
            iCor = ll + 1;
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                //Update Y2[iCor]
                temp = 0;
                for(int i = 0; i<M; i++){
                    temp = temp + Svec[i]*F1[j][ll][i];
                }
                Y2[j][iCor] = Y2[j][iCor] + h*(F1[j][iCor][M-1]-F1[j][ll][M-2]) + h * temp;                        
            }          
            
        }
        //Update F1[0][i]
        for(int i = 0; i<M-1; i++){
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                F1[j][0][i] = F1[j][0][i+1];
            }
        }    

        //Update F1[0][M-1]
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            Y_update[j] = Y2[j][0];        
        }  
        for(int j = 0; j<N_posi;j++){
            F1[2*j][0][M-1] = fx_posi(j,Y_update);
            F1[2*j+1][0][M-1] = fv_posi(j,Y_update);
        }
        for(int j = 0; j<N_nega;j++){
            F1[2*N_posi+2*j][0][M-1] = fx_nega(j,Y_update);
            F1[2*N_posi+2*j+1][0][M-1] = fv_nega(j,Y_update);
        } 

        for(int ll = 0; ll<iStep; ll++){
            iCor = ll + 1;
            for(int i = 0; i<M-1; i++){
                for(int j = 0; j<2*(N_posi+N_nega); j++){
                    F1[j][iCor][i] = F1[j][iCor][i+1]; 
                }                
            }
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                Y_update[j] = Y2[j][iCor];        
            }  
            for(int j = 0; j<N_posi;j++){
                F1[2*j][iCor][M-1] = fx_posi(j,Y_update);
                F1[2*j+1][iCor][M-1] = fv_posi(j,Y_update);
            }
            for(int j = 0; j<N_nega;j++){
                F1[2*N_posi+2*j][iCor][M-1] = fx_nega(j,Y_update);
                F1[2*N_posi+2*j+1][iCor][M-1] = fv_nega(j,Y_update);
            }
        }
        
    }   


    
    tend = MPI_Wtime();
    runtime = tend - tstart; 

    if(rank==0){        
        cout<<"t: "<<runtime << endl;
    }

   // ================== MAIN LOOP FOR TIME ==================
    for(int iTime = 2*M-2; iTime<N+M-1; iTime++){
        if(rank==0){
            //prediction loop
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                temp = 0;
                for(int i = 0; i<M; i++){
                    temp = temp + beta_coeff[i]*F1[j][0][i];
                }            
                Y2[j][0] = Y2[j][0] + h*temp;                   
            }         
        }


        //correction loops up to the second last one
        for(int ll = 0; ll<M-2; ll++){
            iCor = ll + 1;
            if(rank==iCor){
                for(int j = 0; j<2*(N_posi+N_nega); j++){
                    temp = 0;
                    for(int i = 0; i<M; i++){
                        temp = temp + Svec[i]*F1[j][ll][i];
                    }
                    Y2[j][iCor] = Y2[j][iCor] + h*(F1[j][iCor][M-1]-F1[j][ll][M-2]) + h * temp;
                }
            }
        }
        // last correction loop
        if(rank==M-1){
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                temp = 0;
                for(int i = 0; i<M; i++){
                    temp = temp + Svec[i]*F1[j][M-2][i];
                }
                Y2[j][M-1] = Y2[j][M-1] + h*(F2[j]-F1[j][M-2][M-2]) + h*temp;    
            }        
        }

        //~~~~~~~~~~~~~~ Updating Stencil ~~~~~~~~~~~
        //---> updating correctors stencil ~~~~~~~~~~
        if(rank==0){           

            for(int j = 0; j<2*(N_posi+N_nega); j++){
                for(int i = 0; i<M-1; i++){
                    F1[j][0][i] = F1[j][0][i+1];
                }
            }
            
            //Update F1[0][M-1]
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                Y_update[j] = Y2[j][0];        
            }  
            for(int j = 0; j<N_posi;j++){
                F1[2*j][0][M-1] = fx_posi(j,Y_update);
                F1[2*j+1][0][M-1] = fv_posi(j,Y_update);
            }
            for(int j = 0; j<N_nega;j++){
                F1[2*N_posi+2*j][0][M-1] = fx_nega(j,Y_update);
                F1[2*N_posi+2*j+1][0][M-1] = fv_nega(j,Y_update);
            } 
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                MPI_Send(&F1[j][0][M-1],1,MPI_DOUBLE,1, j, MPI_COMM_WORLD);  
                MPI_Recv(&yy[j],1,MPI_DOUBLE, M-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            }
        }
        for(int ll = 1; ll<M-1; ll++){
            if(rank==ll){
                for(int i = 0; i<M-1; i++){
                    for(int j = 0; j<2*(N_posi+N_nega); j++){
                        F1[j][ll][i] = F1[j][ll][i+1];
                        F1[j][ll-1][i] = F1[j][ll-1][i+1];                        
                    }
                }
                
                //Update F1[ll][M-1]
                for(int j = 0; j<2*(N_posi+N_nega); j++){
                    Y_update[j] = Y2[j][ll];        
                }  
                for(int j = 0; j<N_posi;j++){
                    F1[2*j][ll][M-1] = fx_posi(j,Y_update);
                    F1[2*j+1][ll][M-1] = fv_posi(j,Y_update);
                }
                for(int j = 0; j<N_nega;j++){
                    F1[2*N_posi+2*j][ll][M-1] = fx_nega(j,Y_update);
                    F1[2*N_posi+2*j+1][ll][M-1] = fv_nega(j,Y_update);
                }
                for(int j = 0; j<2*(N_posi+N_nega); j++){
                    MPI_Send(&F1[j][ll][M-1],1,MPI_DOUBLE,ll+1, j, MPI_COMM_WORLD);  
                    MPI_Recv(&F1[j][ll-1][M-1],1,MPI_DOUBLE,ll-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }

            }

        }
        // storing the final answer
        if(rank==M-1){
            //Update F2
            for(int j = 0; j<2*(N_posi+N_nega); j++){
                Y_update[j] = Y2[j][M-1];        
            }  

            for(int j = 0; j<N_posi;j++){
                F2[2*j] = fx_posi(j,Y_update);
                F2[2*j+1] = fv_posi(j,Y_update);
            }
            for(int j = 0; j<N_nega;j++){
                F2[2*N_posi+2*j] = fx_nega(j,Y_update);
                F2[2*N_posi+2*j+1] = fv_nega(j,Y_update);
            }    

            for(int i = 0; i<M-1; i++){
                for(int j = 0; j<2*(N_posi+N_nega); j++){
                    F1[j][M-2][i] = F1[j][M-2][i+1];
                }
            }            

            for(int j = 0; j<2*(N_posi+N_nega); j++){
                MPI_Send(&Y2[j][M-1],1,MPI_DOUBLE,0, j, MPI_COMM_WORLD);  
                MPI_Recv(&F1[j][M-2][M-1],1,MPI_DOUBLE,M-2, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            }

            // ---> updating predictor stencil
            // ** approach #0:         
        }
        MPI_Barrier;

    }
    tend2 = MPI_Wtime();
    runtime = tend2 - tstart; 



    if(rank==0){        
        cout<<"t: "<<runtime << endl;
        for(int j = 0; j<2*(N_posi+N_nega); j++){
            cout<<yy[j]<<' ';
        }       
    }

    MPI_Finalize();

    if (rank == 0){          
        fstream myfile;
        string name("RIDC_100");
        string txt(".txt");
        name = name + txt;    
        myfile.open(name,fstream::out);
        for (int i = 0; i < 2*(N_posi+N_nega); i++) {
            myfile << fixed << setprecision(12)<<yy[i] << ",";
            myfile << endl; 
        }    
        myfile.close();
    } 

    return 0;
}