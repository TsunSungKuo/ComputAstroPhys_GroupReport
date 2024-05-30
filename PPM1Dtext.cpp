# include <iostream>
# include <stdio.h>
# include <cmath>
#include <algorithm>

using namespace std;

// --- Initial Parameter --------
const int N_In = 128;
const int nghost = 2;
const int N = N_In + 2*nghost;

double cfl = 1.0;
double gamma_val = 5.0/3.0;
double Leng = 1.0;
double U[N][5] = {0}; //Conservative prop.



double end_time = 0.4;
double t = 0.0;
int limiter = 0;
double dx = Leng/N_In; 
double dt ;

// --- Trick ---
struct vec5
{
    double u[5];
};
struct vecLR
{
    double L[N][N][N][5];
    double R[N][N][N][5];
};
struct vec
{
    double U[N][N][N][5];

};
struct matrix
{
    double mat[N][5];
};

// -------------------------------------------------------------------
//define initial condition
// -------------------------------------------------------------------
vec5 InitialCondition( double x )
{
// Sod shock tube
    double d, u, v, w, P, E;
    if( x < 0.5*Leng ){
        d = 1.25e3;  // density
        u = 0.0;  // velocity x
        v = 0.0;  // velocity y
        w = 0.0;  // velocity z
        P = 5.0e2;  // pressure
        // E = P/(gamma_val-1.0) + 0.5*d*( u*u + v*v + w*w );   // energy density
    }
    else{
        d = 1.25e2;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        P = 5.0;
        // E = P/(gamma_val-1.0) + 0.5*d*( u*u + v*v + w*w );
    }
    //  conserved variables [0/1/2/3/4] <--> [density/momentum x/momentum y/momentum z/energy]
    vec5 dd = {0};
    for(int pos=0; pos<N; pos++)
    {
        //dd.u = {d, d*u, d*v, d*w, E};
        
        dd.u[0] = d;
        dd.u[1] = u;
        dd.u[2] = v;
        dd.u[3] = w;
        dd.u[4] = P;
        
    }
    return dd;
}




void BoundaryCondition( double U[N][5] )
{
    // outflow
    for(int pos1=0; pos1<N; pos1++){

        for(int j=0; j<5; j++){

            if(pos1 < nghost ){
                U[pos1][j] = U[nghost][j];
            }
            else if (pos1 > N-nghost)
            {
                U[pos1][j] = U[N-nghost][j];
            

        }
        
        }
    }
}


double ComputeLimitedSlope(double L, double C, double R )
{
//  compute the left and right slopes
    double slope_L;
    double slope_R;
    double slope_LR;
    double slope_limited;

    double slope_LplusR;
    slope_L = C - L;
    slope_R = R - C;
    slope_LplusR = (slope_L+slope_R)/2.;
    slope_LR = slope_L*slope_R;
        if(slope_LR>0.0 && slope_L>0.0 ){
                slope_limited =std::min({abs(slope_LplusR),abs(slope_L),abs(slope_R)});
            }else if(slope_LR>0.0 && slope_L<0.0 ){
                slope_limited = -std::min({abs(slope_LplusR),abs(slope_L),abs(slope_R)});
            }else{
                slope_limited = 0.0;
            }
    return slope_limited;
}

void Computehalf(double a_jl,double a_j, double a_jr, double slopel,double slope,double sloper,double a_jhalf[2] )
{

    double a_jhr;
    double a_jhl;
        a_jhr = a_j+(a_jr-a_j)/2.+1./6.*(slope-sloper);
        a_jhl = a_jl+(a_j-a_jl)/2.+1./6.*(slopel-slope);
    if((a_jhr-a_j)*(a_j-a_jhl)<=0.){
        a_jhr = a_j;
        a_jhl = a_j;
    }
    else if((a_jhr-a_jhl)*(a_j-0.5*(a_jhl+a_jhr))>pow((a_jhr-a_jhl),2)/6.){
        a_jhl = 3.*a_j-2.*a_jhr;
    }
    else if((a_jhr-a_jhl)*(a_j-0.5*(a_jhl+a_jhr))<-pow((a_jhr-a_jhl),2)/6.){
        a_jhr = 3.*a_j-2.*a_jhl;
    }
    a_jhalf[0]=a_jhl;
    a_jhalf[1]=a_jhr;
}
void DataReconstruction_PPM( double U[N][5],double dt,double gamma_val,double dx )
{
    double U_n1[N][5]={0};
    double slope_x[N][5]={0};

    double hvalue_x[N][5][2]={0};

    double d_dxyz[N][5]={0};
    double calmax[2];
    double a_jhalf[2]={0};


 
//  allocate memory

    for(int num=0; num<5; num++){

    for(int k=1; k<N-1; k++){
        slope_x[k][num] = ComputeLimitedSlope( U[k-1][num], U[k][num], U[k+1][num] );
  
    }
    
    for(int k=1; k<N-1; k++){
    Computehalf(U[k-1][num], U[k][num], U[k+1][num], slope_x[k-1][num], slope_x[k][num], slope_x[k+1][num], a_jhalf);
    hvalue_x[k][num][0]=a_jhalf[0];
    hvalue_x[k][num][1]=a_jhalf[1];


    }
    
    }    
    

    for(int k=2; k<N-2; k++){
    for(int num=0; num<5; num++){
  
    if(U[k][1]>0){
    d_dxyz[k][num]=hvalue_x[k-1][num][1]-hvalue_x[k][num][1]-dt*U[k][1]/dx/2.*((2.-2.*U[k][1]*dt/dx)*(hvalue_x[k-1][num][0]-hvalue_x[k][num][0])+(4.-2.*U[k][1]*dt/dx)*(hvalue_x[k-1][num][1]-hvalue_x[k][num][1])-(6.-4.*U[k][1]*dt/dx)*(U[k-1][num]-U[k][num]));
    }else {
    d_dxyz[k][num]=hvalue_x[k][num][0]-hvalue_x[k+1][num][0]-dt*U[k][1]/dx/2.*((2.+2.*U[k][1]*dt/dx)*(hvalue_x[k][num][0]-hvalue_x[k+1][num][0])+(4.+2.*U[k][1]*dt/dx)*(hvalue_x[k][num][1]-hvalue_x[k+1][num][1])+(6.+4.*U[k][1]*dt/dx)*(U[k][num]-U[k+1][num]));
    }
    // else{
    //     d_dxyz[k][num]=(hvalue_x[k-1][num][1]-hvalue_x[k][num][1]+hvalue_x[k][num][0]-hvalue_x[k+1][num][0])/2.;
    // }
    
    }
    
    }

    for(int i=2; i<N-2; i++){

    //denaity
    U_n1[i][0]=U[i][0]+dt/dx*(U[i][1]*d_dxyz[i][0]+U[i][0]*(d_dxyz[i][1]));
    //velosity_x
    U_n1[i][1]=U[i][1]+dt/dx*(U[i][1]*d_dxyz[i][1]+d_dxyz[i][4]/U[i][0]);
   //pressure
    U_n1[i][4]=U[i][4]+dt/dx*(gamma_val*U[i][4]*(d_dxyz[i][1])+U[i][1]*d_dxyz[i][4]);
     }
    for(int i=2; i<N-2; i++){

    //denaity
    U[i][0]=U_n1[i][0];
    //velosity_x
    U[i][1]=U_n1[i][1];
    //velosity_y
   //pressure
    U[i][4]=U_n1[i][4];
     }
    BoundaryCondition( U );

    for(int i=2; i<N-2; i++){
       calmax[0]=pow(pow(U[i][1],2.)+pow(U[i][2],2.)+pow(U[i][3],2.),0.5)+gamma_val*U[i][4]/U[i][0];
    if(calmax[0]>calmax[1]){
        calmax[1]=calmax[0];}
    }
    dt = dx/2./calmax[1];
    

}


int main()
{
    double calmax[2];
    FILE *output;
    output = fopen("output.txt","w");
    double x;
    double U[N][5];
    for(int pos=0; pos<N; pos++)
    {
        vec5 dd = {0};
        x = (pos + 0.5) * dx;
        if(pos>=nghost and pos<=N-nghost)
        {
            dd = InitialCondition(x);
        }
        
        for(int i=0; i<5; i++)
        {
            U[pos][i] = dd.u[i];
        }
    }
    for(int i=2; i<N-2; i++){
       calmax[0]=pow(pow(U[i][1],2.)+pow(U[i][2],2.)+pow(U[i][3],2.),0.5)+gamma_val*U[i][4]/U[i][0];
    if(calmax[0]>calmax[1]){
        calmax[1]=calmax[0];}
    }
    dt = dx/2./calmax[1];
    while (t >=0 )
    {



        DataReconstruction_PPM( U, dt, gamma_val, dx );
        t=t+dt;
        fprintf( output, "t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);
        printf("t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);
    
    for(int pos1=0; pos1<N; pos1++){

            // THIS LOOP IS FOR PRINT TO DATA
            double data_conservative[5] = {U[pos1][0],U[pos1][1],U[pos1][2],U[pos1][3],U[pos1][4]};
            
            
            for(int i=0; i<5; i++)
                fprintf(output, "%.9e ",data_conservative[i]);
            fprintf(output, "\n");
        }

        // update time
        if(t >= end_time)
            break;
    }
    fclose(output);
}
