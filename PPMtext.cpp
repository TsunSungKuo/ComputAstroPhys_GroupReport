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
double gamma = 5.0/3.0;
double Leng = 1.0;
double U[N][5] = {0}; //Conservative prop.



double end_time = 0.4;
double t = 0.0;
int limiter = 0;
double dx = Leng/N_In; 
double dt;

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
        // E = P/(gamma-1.0) + 0.5*d*( u*u + v*v + w*w );   // energy density
    }
    else{
        d = 1.25e2;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        P = 5.0;
        // E = P/(gamma-1.0) + 0.5*d*( u*u + v*v + w*w );
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


// vec5 Conserved2Primitive( double U[5] )
// {
//     vec5 W = {0};

//     W.u[0] = U[0];
//     W.u[1] = U[1]/U[0];
//     W.u[2] = U[2]/U[0];
//     W.u[3] = U[3]/U[0];
//     W.u[4] = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
    
//     //printf("c2p = %f,%f,%f,%f,%f\n",W.u[0],W.u[1],W.u[2],W.u[3],W.u[4]);

//     return W ;
// }
// -------------------------------------------------------------------
// convert primitive variables to conserved variables
// -------------------------------------------------------------------
// vec5 Primitive2Conserved( double W[5] )
// {
//     vec5 U = {0};

//     U.u[0] = W[0];
//     U.u[1] = W[0]*W[1];
//     U.u[2] = W[0]*W[2];
//     U.u[3] = W[0]*W[3];
//     U.u[4] = W[4]/(gamma_val-1.0) + 0.5*W[0]*( W[1]*W[1] + W[2]*W[2] + W[3]*W[3] );

//     //printf("p2c = %f,%f,%f,%f,%f\n",U.u[0],U.u[1],U.u[2],U.u[3],U.u[4]);
//     return U;
// }

// vec5 Conserved2Flux( double U[5] )
// {
//     vec5 flux = {0};

//     double P = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
    
//     double u = {1};
//     u = U[1] / U[0];

//     flux.u[0] = U[1];
//     flux.u[1] = u*U[1] + P;
//     flux.u[2] = u*U[2];
//     flux.u[3] = u*U[3];
//     flux.u[4] = u*( U[4] + P );


//     //printf("C2F = %f,%f,%f,%f,%f\n",flux.u[0],flux.u[1],flux.u[2],flux.u[3],flux.u[4]);
//     return flux;
// }

void BoundaryCondition( double U[N][N][N][5] )
{
    // outflow
    for(int pos1=0; pos1<N; pos1++){
        for(int pos2=0; pos2<N; pos2++){
            for(int pos3=0; pos3<N; pos3++){
        for(int j=0; j<5; j++){

            if(pos1 < nghost ){
                U[pos1][pos2][pos3][j] = U[nghost][pos2][pos3][j];
            }
            else if (pos1 > N-nghost)
            {
                U[pos1][pos2][pos3][j] = U[N-nghost][pos2][pos3][j];
            
            } else if (pos2 > N-nghost)
            {
                U[pos1][pos2][pos3][j] = U[pos1][N-nghost][pos3][j];
            
            }else if (pos3 > N-nghost)
            {
                U[pos1][pos2][pos3][j] = U[pos1][pos2][N-nghost][j];
            
            }
            else if (pos3 > N-nghost)
            {
                U[pos1][pos2][pos3][j] = U[pos1][pos2][N-nghost][j];
            
            }
        }
        }
        }
    }
}

// -------------------------------------------------------------------
// compute time-step by the CFL condition
// -------------------------------------------------------------------
// double ComputeTimestep( double U[N][5] )
// {
//     double P[N], a[N], u[N], v[N], w[N];
//     double ua_max=0;
//     for(int pos=0; pos<N; pos++)
//     {
//         P[pos] = ComputePressure( U[pos][0], U[pos][1], U[pos][2], U[pos][3], U[pos][4] );
//         a[pos] = sqrt( gamma_val*P[pos]/U[pos][0] );
//         u[pos] = abs( U[pos][1]/U[pos][0] );
//         v[pos] = abs( U[pos][2]/U[pos][0] );
//         w[pos] = abs( U[pos][3]/U[pos][0] );

//         if(u[pos]+a[pos]>ua_max)    //find maximum
//             ua_max = u[pos] + a[pos];
//     }
//     double max_info_speed = ua_max;
//     //max_info_speed = std::max( ua )
//     double dt_cfl         = cfl*dx/max_info_speed;
//     double dt_end         = end_time - t;

//     return min( dt_cfl, dt_end );
// }

// -------------------------------------------------------------------
// compute limited slope
// -------------------------------------------------------------------

// vec5 Roe( double L[5], double R[5] )
// {
//     vec5 flux;

//     //compute the enthalpy of the left and right states: H = (E+P)/rho
//     double P_L = ComputePressure( L[0], L[1], L[2], L[3], L[4] );
//     double P_R = ComputePressure( R[0], R[1], R[2], R[3], R[4] );
//     double H_L = ( L[4] + P_L )/L[0];
//     double H_R = ( R[4] + P_R )/R[0];

//     //compute Roe average values
//     double rhoL_sqrt = sqrt(L[0]);
//     double rhoR_sqrt = sqrt(R[0]);

//     double u  = ( L[1]/rhoL_sqrt + R[1]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
//     double v  = ( L[2]/rhoL_sqrt + R[2]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
//     double w  = ( L[3]/rhoL_sqrt + R[3]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
//     double H  = ( rhoL_sqrt*H_L  + rhoR_sqrt*H_R  ) / ( rhoL_sqrt + rhoR_sqrt );
//     double V2 = u*u + v*v + w*w;
//     //check negative pressure
//     //assert H-0.5*V2 > 0.0, "negative pressure!"
//     double a  = sqrt( (gamma_val-1.0)*(H - 0.5*V2));

//     //compute the amplitudes of different characteristic waves
//     double dU[5];
//     for(int i=0; i<5; i++){
//         dU[i] = R[i] - L[i];
//     }
//     double amp[5] = {0};
//     amp[2] = dU[2] - v*dU[0];
//     amp[3] = dU[3] - w*dU[0];
//     amp[1] = (gamma_val-1.0)/(a*a)*( dU[0]*(H-u*u) + u*dU[1] - dU[4] + v*amp[2] + w*amp[3] );
//     amp[0] = (0.5/a)*( dU[0]*(u+a) - dU[1] - a*amp[1] );
//     amp[4] = dU[0] - amp[0] - amp[1];


//     // compute the eigenvalues and right eigenvector matrix
//     double EigenValue[5]    = {u-a, u, u, u, u+a};
//     double EigenVector_R[5][5] = {{1.0, u-a,   v,   w,  H-u*a},
//                                   {1.0,   u,   v,   w, 0.5*V2},
//                                   {0.0, 0.0, 1.0, 0.0,      v},
//                                   {0.0, 0.0, 0.0, 1.0,      w},
//                                   {1.0, u+a,   v,   w,  H+u*a}};

//     //compute the fluxes of the left and right states
//     vec5 flux_L_temp = Conserved2Flux( L );
//     vec5 flux_R_temp = Conserved2Flux( R );
    
//     double flux_L[5] = {0};
//     double flux_R[5] = {0};
//     for(int i=0;i<5;i++)
//     {
//         flux_L[i] = flux_L_temp.u[i];
//         flux_R[i] = flux_R_temp.u[i];
//     }
    
//     //compute the Roe flux
//     for(int i=0; i<5; i++){
//         amp[i] = amp[i]*abs(EigenValue[i]);
//     } 

//     double temp[5] ={0}; //a temporary storeing matrix multipycation 
//     for(int i=0; i<5; i++){
//         for(int j=0; j<5; j++){
//             temp[i] += amp[j]*EigenVector_R[j][i];
//         }
//     }
    
//     for(int i=0; i<5; i++)
//     {
//         flux.u[i] = 0.5*( flux_L[i] + flux_R[i] ) - 0.5*temp[i];
//     }

//     return flux;
    
// }
double ComputeLimitedSlope(double L, double C, double R )
{
//  compute the left and right slopes
    double slope_L;
    double slope_R;
    double slope_LR;
    double slope_limited;


    slope_L = C - L;
    slope_R = R - C;
    slope_LR = slope_L*slope_R;
        if(slope_LR>0.0 && slope_L>0.0 ){
                slope_limited =min(abs(slope_L+slope_R)/2.,abs(slope_L),abs(slope_R));
            }else if(slope_LR>0.0 && slope_L<0.0 ){
                slope_limited = 0-min(abs(slope_L+slope_R)/2.,abs(slope_L),abs(slope_R));
            }else{
                slope_limited = 0.0;
            }
    return slope_limited;
}

void Computehalf(double a_jl,double a_j, double a_jr, double slopel,double slope,double sloper,double a_jhalf[2] )
{

    double a_jhr;
    double a_jhl;
        a_jhr = a_j+(a_jr-a_j)/2.+1/6.*(slope-sloper);
        a_jhl = a_jl+(a_j-a_jl)/2.+1/6.*(slopel-slope);
    if((a_jhr-a_j)*(a_j-a_jhl)<=0.){
        a_jhr = a_j;
        a_jhl = a_j;
    }
    else if((a_jhr-a_j)*(a_j-0.5*(a_jhl+a_jhr))>pow((a_jhr-a_jhl),2)/6.){
        a_jhl = 3.*a_j-2.*a_jhr;
    }
    else if((a_jhr-a_j)*(a_j-0.5*(a_jhl+a_jhr))<pow((a_jhr-a_jhl),2)/6.){
        a_jhr = 3.*a_j-2.*a_jhl;
    }
    a_jhalf[0]=a_jhl;
    a_jhalf[1]=a_jhr;
}
void DataReconstruction_PPM( double U[N][N][N][5],double dt,double gamma,double dx )
{
//  allocate memory
    double U_n1[N][N][N][5]={0};
    double slope_x[N][N][N][5]={0};
    double slope_y[N][N][N][5]={0};
    double slope_z[N][N][N][5]={0};
    double hvalue_x[N][N][N][5][2]={0};
    double hvalue_y[N][N][N][5][2]={0};
    double hvalue_z[N][N][N][5][2]={0};
    double d_dxyz[N][N][N][5][3]={0};
    double calmax[2];
    double a_jhalf[2]={0};
    for(int num=0; num<5; num++){
    for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
    for(int k=1; k<N-1; k++){
        slope_x[k][i][j][num] = ComputeLimitedSlope( U[k-1][i][j][num], U[k][i][j][num], U[k+1][i][j][num] );
        slope_y[j][k][i][num] = ComputeLimitedSlope( U[j][k-1][i][num], U[j][k][i][num], U[j][k+1][i][num] );
        slope_z[i][j][k][num] = ComputeLimitedSlope( U[i][j][k-1][num], U[i][j][k][num], U[i][j][k+1][num] );
    }
    }
    }
    for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
    for(int k=1; k<N-1; k++){
    Computehalf(U[k-1][i][j][num], U[k][i][j][num], U[k+1][i][j][num], slope_x[k-1][i][j][num], slope_x[k][i][j][num], slope_x[k+1][i][j][num], a_jhalf);
    hvalue_x[k][i][j][num][0]=a_jhalf[0];
    hvalue_x[k][i][j][num][1]=a_jhalf[1];
    Computehalf(U[j][k-1][i][num], U[j][k][i][num], U[j][k+1][i][num], slope_y[j][k-1][i][num], slope_y[j][k][i][num], slope_y[j][k+1][i][num], a_jhalf);
    hvalue_y[j][k][i][num][0]=a_jhalf[0];
    hvalue_y[j][k][i][num][1]=a_jhalf[1];
    Computehalf(U[i][j][k-1][num], U[i][j][k][num], U[i][j][k+1][num], slope_z[i][j][k-1][num], slope_z[i][j][k][num], slope_z[i][j][k+1][num], a_jhalf);
    hvalue_z[i][j][k][num][0]=a_jhalf[0];
    hvalue_z[i][j][k][num][1]=a_jhalf[1];
    }
    }
    }    
    }
    for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
    for(int k=2; k<N-2; k++){
    for(int num=0; num<5; num++){
    for(int num1=0; num1<3; num1++){
    if(U[k][i][j][1]>0){
    d_dxyz[k][i][j][num][0]=hvalue_x[k-1][i][j][num][1]-hvalue_x[k][i][j][num][1]-dt*U[k][i][j][1]/2.*((2.-2.*U[k][i][j][1]*dt)*(hvalue_x[k-1][i][j][num][0]-hvalue_x[k][i][j][num][0])+(4.-2.*U[k][i][j][1]*dt)*(hvalue_x[k-1][i][j][num][1]-hvalue_x[k][i][j][num][1])-(6.-4.*U[k][i][j][1]*dt)*(U[k-1][i][j][num]-U[k][i][j][num]));
    }else{
    d_dxyz[k][i][j][num][0]=hvalue_x[k][i][j][num][1]-hvalue_x[k+1][i][j][num][1]-dt*U[k][i][j][1]/2.*((2.-2.*U[k][i][j][1]*dt)*(hvalue_x[k][i][j][num][0]-hvalue_x[k+1][i][j][num][0])+(4.-2.*U[k][i][j][1]*dt)*(hvalue_x[k][i][j][num][1]-hvalue_x[k+1][i][j][num][1])+(6.+4.*U[k][i][j][1]*dt)*(U[k][i][j][num]-U[k+1][i][j][num]));
    }
    if(U[j][k][i][2]>0){
    d_dxyz[j][k][i][num][1]=hvalue_y[j][k-1][i][num][1]-hvalue_y[j][k][i][num][1]-dt*U[j][k][i][1]/2.*((2.-2.*U[j][k][i][1]*dt)*(hvalue_y[j][k-1][i][num][0]-hvalue_y[j][k][i][num][0])+(4.-2.*U[j][k][i][1]*dt)*(hvalue_y[j][k-1][i][num][1]-hvalue_y[j][k][i][num][1])-(6.-4.*U[j][k][i][1]*dt)*(U[j][k-1][i][num]-U[j][k][i][num]));
    }else{
    d_dxyz[j][k][i][num][1]=hvalue_y[j][k][i][num][1]-hvalue_y[j][k+1][i][num][1]-dt*U[j][k][i][1]/2.*((2.-2.*U[j][k][i][1]*dt)*(hvalue_y[j][k][i][num][0]-hvalue_y[j][k+1][i][num][0])+(4.-2.*U[j][k][i][1]*dt)*(hvalue_y[j][k][i][num][1]-hvalue_y[j][k+1][i][num][1])+(6.+4.*U[j][k][i][1]*dt)*(U[j][k][i][num]-U[j][k+1][i][num]));
    }
    if(U[i][j][k][3]>0){
    d_dxyz[i][j][k][num][2]=hvalue_z[i][j][k-1][num][1]-hvalue_z[i][j][k][num][1]-dt*U[i][j][k][1]/2.*((2.-2.*U[i][j][k][1]*dt)*(hvalue_z[i][j][k-1][num][0]-hvalue_z[i][j][k][num][0])+(4.-2.*U[i][j][k][1]*dt)*(hvalue_z[i][j][k-1][num][1]-hvalue_z[i][j][k][num][1])-(6.-4.*U[i][j][k][1]*dt)*(U[i][j][k-1][num]-U[i][j][k][num]));
    }else{
    d_dxyz[i][j][k][num][2]=hvalue_z[i][j][k][num][1]-hvalue_z[i][j][k+1][num][1]-dt*U[i][j][k][1]/2.*((2.-2.*U[i][j][k][1]*dt)*(hvalue_z[i][j][k][num][0]-hvalue_z[i][j][k+1][num][0])+(4.-2.*U[i][j][k][1]*dt)*(hvalue_z[i][j][k][num][1]-hvalue_z[i][j][k+1][num][1])+(6.+4.*U[i][j][k][1]*dt)*(U[i][j][k][num]-U[i][j][k+1][num]));
    }
    }
    }
    }
    }
    }
    for(int i=2; i<N-2; i++){
    for(int j=2; j<N-2; j++){
    for(int k=2; k<N-2; k++){

    //denaity
    U_n1[i][j][k][0]=U[i][j][k][0]-dt/dx*(U[i][j][k][1]*d_dxyz[i][j][k][0][0]+U[i][j][k][2]*d_dxyz[i][j][k][0][1]+U[i][j][k][2]*d_dxyz[i][j][k][0][2]+U[i][j][k][0]*(d_dxyz[i][j][k][1][0]+d_dxyz[i][j][k][2][1]+d_dxyz[i][j][k][3][2]));
    //velosity_x
    U_n1[i][j][k][1]=U[i][j][k][1]-dt/dx*(U[i][j][k][1]*d_dxyz[i][j][k][1][0]+U[i][j][k][2]*d_dxyz[i][j][k][1][1]+U[i][j][k][3]*d_dxyz[i][j][k][1][2]+d_dxyz[i][j][k][4][0]/U[i][j][k][0]);
    //velosity_y
    U_n1[i][j][k][2]=U[i][j][k][2]-dt/dx*(U[i][j][k][1]*d_dxyz[i][j][k][2][0]+U[i][j][k][2]*d_dxyz[i][j][k][2][1]+U[i][j][k][3]*d_dxyz[i][j][k][2][2]+d_dxyz[i][j][k][4][1]/U[i][j][k][0]);
    //velosity_z
    U_n1[i][j][k][3]=U[i][j][k][3]-dt/dx*(U[i][j][k][1]*d_dxyz[i][j][k][3][0]+U[i][j][k][2]*d_dxyz[i][j][k][3][1]+U[i][j][k][3]*d_dxyz[i][j][k][3][2]+d_dxyz[i][j][k][4][2]/U[i][j][k][0]);
    //pressure
    U_n1[i][j][k][4]=U[i][j][k][4]-dt/dx*(gamma*U[i][j][k][4]*(d_dxyz[i][j][k][1][0]+d_dxyz[i][j][k][2][1]+d_dxyz[i][j][k][3][2])+U[i][j][k][1]*d_dxyz[i][j][k][4][0]+U[i][j][k][2]*d_dxyz[i][j][k][4][1]+U[i][j][k][3]*d_dxyz[i][j][k][4][2]);
    }
    }
    }
    BoundaryCondition( U );
    for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
    for(int k=0; k<N; k++){

    calmax[0]=pow(pow(U[i][j][k][1],2.)+pow(U[i][j][k][2],2.)+pow(U[i][j][k][3],2.),0.5)+gamma*U[i][j][k][4];
    if(calmax[0]>calmax[1]){calmax[1]=calmax[0];}
    }
    }
    }
    dt= 1./calmax[1];

}

//----------------------------------------------------------------
/*
int main()
{
    double a[5] = {1.25e3, 0.0, 0.0, 0.0, 5.0e2};
    vec5 b = Primitive2Conserved(a);
    vec5 c = Conserved2Primitive(b.u);
    for(int i=0; i<5; i++)
    printf("%f  ",c.u[i]); 
}
*/
// double Calculate_slope(double value[N][N][N] )
// {
//     double W[N][N][N][3]={0};
//     //slope in x-axis
//     for(int i=1; i<N-1; i++)
//     for(int j=0; j<N; j++)
//     for(int k=0; k<N; k++)
//     double W[i][j][k][0]= std::min((value[N+1][N][N]-value[N-1][N][N])/2 ,(value[N+1][N][N]-value[N][N][N]),(value[N][N][N]-value[N-1][N][N])/2) ;
//     ;;;
//     //slope in y-axis
//     for(int i=1; i<N-1; i++)
//     for(int j=0; j<N; j++)
//     for(int k=0; k<N; k++)
//     double W[i][j][k][1]= std::min((value[N][N+1][N]-value[N][N-1][N])/2 ,(value[N][N+1][N]-value[N][N][N]),(value[N][N][N]-value[N][N-1][N])/2) ;
//     ;;;
//     //slope in z-axis
//     for(int i=1; i<N-1; i++)
//     for(int j=0; j<N; j++)
//     for(int k=0; k<N; k++)
//     double W[i][j][k][2]= std::min((value[N][N][N+1]-value[N][N][N-1])/2 ,(value[N][N][N+1]-value[N][N][N]),(value[N][N][N]-value[N][N][N-1])/2) ;
//     ;;;
     
//     //assert np.all( P > 0 ), "negative pressure !!"
    
//     return W;
//}
//double Calculate_half(double value[N][N][N] )

int main()
{
    FILE *output;
    output = fopen("output.txt","w");
    double x;
    double U[N][N][N][5];
   for(int pos1=0; pos1<N; pos1++)
    for(int pos2=0; pos2<N; pos1++)
    for(int pos3=0; pos3<N; pos1++)
    {
        {
            {
        vec5 dd = {0};
        x = (pos1 + 0.5) * dx;
        if((pos1>=nghost and pos1<=N-nghost )||( pos2>=nghost and pos2<=N-nghost  )||( pos2>=nghost and pos2<=N-nghost))
        {
            dd = InitialCondition(x);
        }

        for(int i=0; i<5; i++)
        {
            U[pos1][pos2][pos3][i] = dd.u[i];
        }
    }
    }
    }
    
    while (t >=0 )
    {
        // Set boundary condition


        // calculate time-step



        DataReconstruction_PPM( U, dt, gamma, dx );
        t=t+dt;
        fprintf( output, "t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);
        printf("t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);
        //updat the face-center variales by 0.5*dt
        // for(int pos=1; pos<N-1; pos++)
        // {
        //     double L_temp[5], R_temp[5];
        //     for(int i=0; i<5; i++)
        //     {
        //         L_temp[i] = LR.L[pos][i];
        //         R_temp[i] = LR.R[pos][i];
        //         //printf("LR.L=%f\n",LR.L[pos][i]);
        //     }

        //     vec5 Flux_L, Flux_R;
        //     Flux_L = Conserved2Flux(L_temp);
        //     Flux_R = Conserved2Flux(R_temp);
            
        //     double dflux[5] = {0};
        //     for(int i=0; i<5; i++)
        //     {
        //     dflux[i] = 0.5*dt/dx*(Flux_R.u[i]-Flux_L.u[i]);
        //     //printf("dflux[%d][%d]=%f\n",pos,i,dflux[i]);
        //     LR.L[pos][i] -= dflux[i];
        //     LR.R[pos][i] -= dflux[i];
        //     }
        // }

        //printf("debug : compute fluxes\n");
        //compute fluxes
        // double flux[N][5];
        // for(int pos=nghost; pos<N-nghost+1; pos++)
        // {
        //     double R_temp[5], L_temp[5];
        //     for(int i=0; i<5; i++)
        //     {
        //         R_temp[i] = LR.R[pos-1][i];
        //         L_temp[i] = LR.L[pos][i];
        //     }
        //     // R_temp is the LEFT state at the j+1/2 interface
        //     vec5 flux_temp = {0};
        //     flux_temp = Roe(R_temp, L_temp);
        //     for(int i=0; i<5; i++)
        //     {
        //         flux[pos][i] = flux_temp.u[i];
        //         //printf("flux[%d][%d]=%f\n",pos,i,flux[pos][i]);
        //     }
        // }

        // update the volume-averaged input variables by dt
    for(int pos1=0; pos1<N; pos1++){
        for(int pos2=0; pos2<N; pos2++){
            for(int pos3=0; pos3<N; pos3++){
            // THIS LOOP IS FOR PRINT TO DATA
            double data_conservative[5] = {U[pos1][pos2][pos3][0],U[pos1][pos2][pos3][1],U[pos1][pos2][pos3][2],U[pos1][pos2][pos3][3],U[pos1][pos2][pos3][4]};
            
            
            for(int i=0; i<5; i++)
                fprintf(output, "%.9e ",data_conservative[i]);
            fprintf(output, "\n");
        }
  
    }

    }
        // update time
        if(t >= end_time)
            break;
    }
    fclose(output);
}
