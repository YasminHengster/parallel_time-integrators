#include <iostream>
#include <math.h>
#include <mpi.h>
#include <iomanip>      // std::setprecision
using namespace std;

double func(double t, double y){
 double f = 4*t*sqrt(y);
 return f;
}


int main(){
  int rank;
  int M=4;

  double T=96; //end time
  double y0=1; //inital condition
  int N=1000; //number of points in the interval [0,T]
  double h = T/N; //stepsize
  int Mm = M-1; //no of corrections
  double yy[N+1]; //array to solution
  yy[0]=y0;
  double t[N+1]; //time array
  double t_ext[N+1+M]; //extended time array
  double F1[Mm][M];
  double F2;
  double Y2[M];
  double Svec[M];
//M=2
//  double S[Mm][Mm+1]={{0.5, 0.5}};
//  double beta[M]= {-1./2, 3./2};
//M=3
//  double S[Mm][Mm+1]={{ 0.41666667,  0.66666667, -0.08333333}, {-0.08333333,  0.66666667,  0.41666667}};
//  double beta[M]={5./12, -16./12, 23./12};
//M=4 
  double S[Mm][Mm+1]={{ 0.375,       0.79166667, -0.20833333,  0.04166667},{-0.04166667,  0.54166667,  0.54166667, -0.04166667}, {0.04166667, -0.20833333,  0.79166667,  0.375     }};
  double beta[M]={-9./24, 37./24, -59./24, 55./24};
//M=5
//  double S[Mm][Mm+1]={{ 0.34861111,  0.89722222, -0.36666667,  0.14722222, -0.02638889},{-0.02638889,  0.48055556,  0.63333333, -0.10277778,  0.01527778},{0.01527778, -0.10277778,  0.63333333,  0.48055556, -0.02638889},{-0.02638889,  0.14722222, -0.36666667, 0.89722222,  0.34861111}};
//  double beta[M]= {251./720, -1274./720, 2616./720, -2774./720, 1901./720};
//M=6
//  double S[Mm][Mm+1]={{ 0.32986111,  0.99097222, -0.55416667,  0.33472222, -0.12013889,  0.01875}, {-0.01875, 0.44236111,  0.70972222, -0.17916667,  0.05347222, -0.00763889},{0.00763889, -0.06458333,  0.55694444,  0.55694444, -0.06458333,  0.00763889},{-0.00763889,  0.05347222, -0.17916667,  0.70972222,  0.44236111, -0.01875},{0.01875,    -0.12013889,  0.33472222, -0.55416667,  0.99097222,  0.32986111}};
//  double beta[M]={-475./720, 2877./720, -7298./720, 9982./720, -7923./720, 4277./720};
  
  //initialise MPI

  MPI_Init(NULL,NULL);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &M);

  //calculate the initial parts in process 0
  if (rank==0){
    //svec:
    for (int i=0; i<M; i++){
    Svec[i]=S[Mm-1][i];
    }	    

    // Value of RHS at initial time
    double F0 = func(0,y0);

    // the time vector
    // float t[N+1];
    for(int i = 0; i<N+1; i++){
        t[i] = i*h;
    }
    std::cout <<t[0] << ";" << std::setprecision(15) << yy[0] << '\n';
    //cout << t[0] << ";" << yy[0] << endl; 
    // extended time vector
    // float t_ext[N+M+1];
    for(int i = 0; i<N+M+1; i++){
        t_ext[i] = i*h;
    }
   
    // float F1[Mm][M];
    for(int i = 0; i<Mm; i++){
        F1[i][0] = F0;
    }
    double F2 = F0;

    //Y2 [M] new point derived in each level
   // float Y2[M];
    for(int i = 0; i<M; i++){
        Y2[i] = y0;
    }	
	
	/* ================== INITIAL PART (1) ==================
    # for this part the predictor and correctors step up to M points in time
    # ** predictor ** uses Runge-Kutta 4 */
    double KK1, KK2, KK3, KK4;
    for(int iTime = 0; iTime<M-1; iTime++){
        KK1 = F1[0][iTime];
        KK2 = func(t[iTime]+h/2, Y2[0]+KK1*h/2);
        KK3 = func(t[iTime]+h/2, Y2[0]+KK2*h/2);
        KK4 = func(t[iTime]+h,   Y2[0]+KK3*h);
        Y2[0] = Y2[0] + h*(KK1 + 2*KK2 + 2*KK3 + KK4)/6;
        F1[0][iTime+1] = func(t[iTime+1], Y2[0]);
    }

    // ** correctors ** use Integral Deffered Correction
    for(int iCor = 1; iCor<M-1; iCor++){
        int ll = iCor - 1;
        for(int iTime = 0; iTime<M-1; iTime++){
            double temp = 0;
            for(int i = 0; i<M; i++){
                temp = temp + S[iTime][i]*F1[ll][i];
            }
            Y2[iCor] = Y2[iCor]+ h*(F1[iCor][iTime]-F1[ll][iTime]) + h * temp;
            F1[iCor][iTime+1] = func(t[iTime+1], Y2[iCor]);
        }
    }

    // treat the last correction loop a little different
    for(int iTime = 0; iTime<M-1; iTime++){
        double temp = 0;
        for(int i = 0; i<M; i++){
            temp = temp + (S[iTime][i]*F1[M-2][i]);
        }
        Y2[M-1] = Y2[M-1] + h*(F2-F1[M-2][iTime]) + h * temp;
        F2 = func(t[iTime+1], Y2[M-1]);
        yy[iTime+1] = Y2[M-1];
    std::cout <<t[iTime+1] << ";" << std::setprecision(15) << yy[iTime+1] << '\n';
	    //cout <<t[iTime+1]<< ";" <<  yy[iTime+1] << endl;
    }


    //================== INITIAL PART (2) ==================
    int iStep, iCor;
    for(int iTime = M-1; iTime<2*M-2; iTime++){
        iStep = iTime - (M-1);

	// prediction loop
        double temp = 0;
        for(int i = 0; i<M; i++){
            temp = temp + beta[i]*F1[0][i];
        }
        Y2[0] = Y2[0] + h*temp;

	//correction loop
        for(int ll = 0; ll<iStep; ll++){
            iCor = ll + 1;
            temp = 0;
            for(int i = 0; i<M; i++){
                temp = temp +S[Mm-1][i]*F1[ll][i];
	    }
            Y2[iCor] = Y2[iCor] + h*(F1[iCor][M-1]-F1[ll][M-2]) + h * temp;
        }
        for(int i = 0; i<M-1; i++){
            F1[0][i] = F1[0][i+1];
        }
        F1[0][M-1] = func(t_ext[iTime+1], Y2[0]);
        for(int ll = 0; ll<iStep; ll++){
            iCor = ll + 1;
            for(int i = 0; i<M-1; i++){
                F1[iCor][i] = F1[iCor][i+1];
	    }
            F1[iCor][M-1] = func(t_ext[iTime+1-iCor], Y2[iCor]);
	    
        }
    }
}


MPI_Bcast(t,N+1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(t_ext,N+1+M, MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(F1,Mm*M,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(Y2,M, MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(Svec,M,MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);
double t_comp1=MPI_Wtime();

//time loop
for (int iTime=2*M-2; iTime< N+M-1; iTime++){
	
	//prediction in rank 0:
	if (rank==0){
		double temp = 0;
            	for(int i = 0; i<M; i++){
                temp = temp + beta[i]*F1[0][i];
            }
	Y2[0]=Y2[0] + h*temp;
	for (int i =1; i<M ; i++){
	MPI_Recv(&Y2[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	}

		
  	//last correction loop M-1 (rank = M-1)
	else if(rank == M-1){
		F2=func(t_ext[iTime-(M-1)], Y2[M-1]);
		double temp =0;
		for(int i=0; i<M; i++){
		temp = temp + S[Mm-1][i]*F1[M-2][i];
		}
	        Y2[M-1] = Y2[M-1] + h * (F2-F1[M-2][M-2]) + h * temp;
		yy[iTime+1-(M-1)]=Y2[M-1];
		std::cout <<t[iTime+1-(M-1)] << ";" << std::setprecision(15) << yy[iTime+1-(M-1)] << '\n';
		//cout << t[iTime+1-(M-1)] << ";" << yy[iTime+1-(M-1)] << endl;
		MPI_Send(&Y2[rank], 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}


	//correction loops 1 to M-2
	else{
	int ll=rank-1;
	double temp=0;
	for (int i=0; i<M; i++){
		temp = temp + Svec[i]*F1[ll][i];
	}
	Y2[rank]=Y2[rank] + h*(F1[rank][M-1]-F1[ll][M-2]) + h*temp;
	MPI_Send(&Y2[rank], 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}


MPI_Barrier(MPI_COMM_WORLD);

//updating stencil
for (int i=0; i<M-1; i++){
	F1[rank][i]=F1[rank][i+1];
}
F1[rank][M-1]=func(t_ext[iTime+1-rank],Y2[rank]);

if (rank==0){
MPI_Send(F1[rank],M,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
}
else if(rank==M-1){
MPI_Recv(F1[rank-1],M,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
else{
MPI_Send(F1[rank],M,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
MPI_Recv(F1[rank-1],M,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
MPI_Barrier(MPI_COMM_WORLD);

//end time loop
}
MPI_Barrier(MPI_COMM_WORLD);
if (rank==0){
double t_comp2= MPI_Wtime();
cout << "# time for everything is " << t_comp2-t_comp1 << endl;
}

  //finalize MPI
  MPI_Finalize();
 return 0; 
 }
