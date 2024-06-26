# include <iostream>
# include <stdio.h>
# include <cmath>
# include <algorithm>
# include <omp.h>
#include <time.h>

using namespace std;

// --- Initial Parameter --------
const int N_In = 128;
const int nghost = 3;
const int N = N_In + 2*nghost;

double cfl = 1.0;
double gamma_val = 5.0/3.0;
double Leng = 1.0;
double U[N][5] = {0}; //Conservative prop.

double dx = Leng/N_In; 
int limiter = 0;

double t = 0.0;
double end_time = 0.4;
double dt;
int time_start, time_end;

int N_thread = 4;

// --- Trick ---
struct vec5
{
    double u[5];
};
struct vecLR
{
    double L[N][5];
    double R[N][5];
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
        E = P/(gamma_val-1.0) + 0.5*d*( u*u + v*v + w*w );   // energy density
    }
    else{
        d = 1.25e2;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        P = 5.0;
        E = P/(gamma_val-1.0) + 0.5*d*( u*u + v*v + w*w );
    }
    //  conserved variables [0/1/2/3/4] <--> [density/momentum x/momentum y/momentum z/energy]
    vec5 dd = {0};
    for(int pos=0; pos<N; pos++)
    {
        //dd.u = {d, d*u, d*v, d*w, E};
        
        dd.u[0] = d;
        dd.u[1] = d*u;
        dd.u[2] = d*v;
        dd.u[3] = d*w;
        dd.u[4] = E;
        
    }
    return dd;
}


double ComputePressure(double d, double px, double py, double pz, double e )
{
    double P = (gamma_val-1.0)*( e - 0.5*(px*px + py*py + pz*pz)/d );
    //assert np.all( P > 0 ), "negative pressure !!"
    
    return P;
}

vec5 Conserved2Primitive( double U[5] )
{
    vec5 W = {0};

    W.u[0] = U[0];
    W.u[1] = U[1]/U[0];
    W.u[2] = U[2]/U[0];
    W.u[3] = U[3]/U[0];
    W.u[4] = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
    
    //printf("c2p = %f,%f,%f,%f,%f\n",W.u[0],W.u[1],W.u[2],W.u[3],W.u[4]);

    return W ;
}
// -------------------------------------------------------------------
// convert primitive variables to conserved variables
// -------------------------------------------------------------------
vec5 Primitive2Conserved( double W[5] )
{
    vec5 U = {0};

    U.u[0] = W[0];
    U.u[1] = W[0]*W[1];
    U.u[2] = W[0]*W[2];
    U.u[3] = W[0]*W[3];
    U.u[4] = W[4]/(gamma_val-1.0) + 0.5*W[0]*( W[1]*W[1] + W[2]*W[2] + W[3]*W[3] );

    //printf("p2c = %f,%f,%f,%f,%f\n",U.u[0],U.u[1],U.u[2],U.u[3],U.u[4]);
    return U;
}

vec5 Conserved2Flux( double U[5] )
{
    vec5 flux = {0};

    double P = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
    
    double u = {1};
    u = U[1] / U[0];

    flux.u[0] = U[1];
    flux.u[1] = u*U[1] + P;
    flux.u[2] = u*U[2];
    flux.u[3] = u*U[3];
    flux.u[4] = u*( U[4] + P );


    //printf("C2F = %f,%f,%f,%f,%f\n",flux.u[0],flux.u[1],flux.u[2],flux.u[3],flux.u[4]);
    return flux;
}

void BoundaryCondition( double U[N][5] )
{
    // outflow
    for(int pos=0; pos<N; pos++){
        for(int j=0; j<5; j++){
            if(pos < nghost){
                U[pos][j] = U[nghost][j];
            }
            else if (pos > N-nghost)
            {
                U[pos][j] = U[N - nghost][j];
            }
        }
    }
}

// -------------------------------------------------------------------
// compute time-step by the CFL condition
// -------------------------------------------------------------------
double ComputeTimestep( double U[N][5] )
{
    double P[N], a[N], u[N], v[N], w[N];
    double ua_max=0;
    for(int pos=0; pos<N; pos++)
    {
        P[pos] = ComputePressure( U[pos][0], U[pos][1], U[pos][2], U[pos][3], U[pos][4] );
        a[pos] = sqrt( gamma_val*P[pos]/U[pos][0] );
        u[pos] = abs( U[pos][1]/U[pos][0] );
        v[pos] = abs( U[pos][2]/U[pos][0] );
        w[pos] = abs( U[pos][3]/U[pos][0] );

        if(u[pos]+a[pos]>ua_max)    //find maximum
            ua_max = u[pos] + a[pos];
    }
    double max_info_speed = ua_max;
    //max_info_speed = std::max( ua )
    double dt_cfl         = cfl*dx/max_info_speed;
    double dt_end         = end_time - t;

    return min( dt_cfl, dt_end );
}

// -------------------------------------------------------------------
// compute limited slope
// -------------------------------------------------------------------
vec5 ComputeLimitedSlope(double L[5], double C[5], double R[5] )
{
//  compute the left and right slopes
    double slope_L[5];
    double slope_R[5];
    double slope_LR[5];
    vec5 slope_limited;

    for(int i=0; i<5; i++)
    {
        slope_L[i] = C[i] - L[i];
        slope_R[i] = R[i] - C[i];
    
        if(limiter==0){
            // apply the van-Leer limiter
            slope_LR[i]      = slope_L[i]*slope_R[i];
            if(slope_LR[i]>0.0){
                slope_limited.u[i] = 2.0*slope_LR[i]/(slope_L[i]+slope_R[i]);
            }
            else{
                slope_limited.u[i] = 0.0;
            }
        }
        else{
            //  apply the min-mod limiter
            slope_LR[i]      = slope_L[i]*slope_R[i];
            if(slope_LR[i]>0.0){
                slope_limited.u[i] = (slope_L[i]/abs(slope_L[i]))*min(abs(slope_L[i]),abs(slope_R[i]));
            }
            else{
                slope_limited.u[i] = 0.0;
            }
        }
    }
    return slope_limited;
}


vec5 Roe( double L[5], double R[5] )
{
    vec5 flux;

    //compute the enthalpy of the left and right states: H = (E+P)/rho
    double P_L = ComputePressure( L[0], L[1], L[2], L[3], L[4] );
    double P_R = ComputePressure( R[0], R[1], R[2], R[3], R[4] );
    double H_L = ( L[4] + P_L )/L[0];
    double H_R = ( R[4] + P_R )/R[0];

    //compute Roe average values
    double rhoL_sqrt = sqrt(L[0]);
    double rhoR_sqrt = sqrt(R[0]);

    double u  = ( L[1]/rhoL_sqrt + R[1]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
    double v  = ( L[2]/rhoL_sqrt + R[2]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
    double w  = ( L[3]/rhoL_sqrt + R[3]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
    double H  = ( rhoL_sqrt*H_L  + rhoR_sqrt*H_R  ) / ( rhoL_sqrt + rhoR_sqrt );
    double V2 = u*u + v*v + w*w;
    //check negative pressure
    //assert H-0.5*V2 > 0.0, "negative pressure!"
    double a  = sqrt( (gamma_val-1.0)*(H - 0.5*V2));

    //compute the amplitudes of different characteristic waves
    double dU[5];
    for(int i=0; i<5; i++){
        dU[i] = R[i] - L[i];
    }
    double amp[5] = {0};
    amp[2] = dU[2] - v*dU[0];
    amp[3] = dU[3] - w*dU[0];
    amp[1] = (gamma_val-1.0)/(a*a)*( dU[0]*(H-u*u) + u*dU[1] - dU[4] + v*amp[2] + w*amp[3] );
    amp[0] = (0.5/a)*( dU[0]*(u+a) - dU[1] - a*amp[1] );
    amp[4] = dU[0] - amp[0] - amp[1];


    // compute the eigenvalues and right eigenvector matrix
    double EigenValue[5]    = {u-a, u, u, u, u+a};
    double EigenVector_R[5][5] = {{1.0, u-a,   v,   w,  H-u*a},
                                  {1.0,   u,   v,   w, 0.5*V2},
                                  {0.0, 0.0, 1.0, 0.0,      v},
                                  {0.0, 0.0, 0.0, 1.0,      w},
                                  {1.0, u+a,   v,   w,  H+u*a}};

    //compute the fluxes of the left and right states
    vec5 flux_L_temp = Conserved2Flux( L );
    vec5 flux_R_temp = Conserved2Flux( R );
    
    double flux_L[5] = {0};
    double flux_R[5] = {0};
    for(int i=0;i<5;i++)
    {
        flux_L[i] = flux_L_temp.u[i];
        flux_R[i] = flux_R_temp.u[i];
    }
    
    //compute the Roe flux
    for(int i=0; i<5; i++){
        amp[i] = amp[i]*abs(EigenValue[i]);
    } 

    double temp[5] ={0}; //a temporary storeing matrix multipycation 
    for(int i=0; i<5; i++){
        for(int j=0; j<5; j++){
            temp[i] += amp[j]*EigenVector_R[j][i];
        }
    }
    
    for(int i=0; i<5; i++)
    {
        flux.u[i] = 0.5*( flux_L[i] + flux_R[i] ) - 0.5*temp[i];
    }

    return flux;
    
}

vecLR DataReconstruction_PLM( double U[N][5] )
{
//  allocate memory
    double W[N][5]={0};
    vecLR LR = {0};

    // conserved variables -> primitive variables
    for(int pos=0; pos<N; pos++)
    {   
        double W_ttemp[5];
        for(int i=0; i<5; i++)
            W_ttemp[i] = U[pos][i];
        
        vec5 W_temp;
        W_temp = Conserved2Primitive(W_ttemp);
        for(int i=0; i<5; i++)
            W[pos][i]  = W_temp.u[i];
    }
    
    for(int pos=1; pos<N-1; pos++)
    {    
        // compute the left and right states of each cell
        double W_minus[5], W_equal[5], W_plus[5];
        for(int i=0; i<5; i++)
        {
            W_minus[i] = W[pos-1][i];
            W_equal[i] = W[pos][i];
            W_plus[i]  = W[pos+1][i];
        }

        vec5 slope_limited;
        slope_limited = ComputeLimitedSlope(W_minus, W_equal, W_plus);

        // get the face-centerd variables
        for(int i=0; i<5; i++)
        {
            LR.L[pos][i] = W[pos][i] - 0.5*slope_limited.u[i];
            LR.R[pos][i] = W[pos][i] + 0.5*slope_limited.u[i];  
            
            // ensure face-centered variable lie between nearby volume-average (~cell-centered)
            LR.L[pos][i] = max( LR.L[pos][i], min(W[pos-1][i], W[pos][i]) );
            LR.L[pos][i] = min( LR.L[pos][i], max(W[pos-1][i], W[pos][i]) );
            LR.R[pos][i] = 2.0*W[pos][i] - LR.L[pos][i];

            LR.R[pos][i] = max( LR.R[pos][i], min(W[pos+1][i], W[pos][i]) );
            LR.R[pos][i] = min( LR.R[pos][i], max(W[pos+1][i], W[pos][i]) );
            LR.L[pos][i] = 2.0*W[pos][i] - LR.R[pos][i];
        }
        // primitve variables -> conserved variables
        
        double L_ttemp[5] = {0};
        double R_ttemp[5] = {1};

        for(int i=0; i<5; i++) 
        {
            L_ttemp[i] = LR.L[pos][i];
            R_ttemp[i] = LR.R[pos][i];
        }

        vec5 L_temp = Primitive2Conserved(L_ttemp);
        vec5 R_temp = Primitive2Conserved(R_ttemp);
       
        for(int i=0; i<5; i++)
        {
            LR.L[pos][i] = L_temp.u[i];
            LR.R[pos][i] = R_temp.u[i];
            //printf("LR.L=%f\n",LR.L[pos][i]);
            //printf("LR.R=%f\n",LR.R[pos][i]); //Promblems before here.
        }
    
    }
    return LR;
}

vecLR DataReconstruction_PPM( double U[N][5] )
{
    //allocate memory
    double W[N][5] = {0};
    vecLR LR = {0};

    // conserved variables -> primitive variables
    for(int pos=0; pos<N; pos++)
    {   
        double W_ttemp[5];
        for(int i=0; i<5; i++)
            W_ttemp[i] = U[pos][i];
        
        vec5 W_temp;
        W_temp = Conserved2Primitive(W_ttemp);
        for(int i=0; i<5; i++)
            W[pos][i]  = W_temp.u[i];
    }

    for(int pos=2; pos<N-2; pos++)
    {
        //compute the left and right states of each cell
        double W_minus[5], W_equal[5], W_plus[5];
        for(int i=0; i<5; i++)
        {
            W_minus[i] = W[pos-1][i];
            W_equal[i] = W[pos][i];
            W_plus[i]  = W[pos+1][i];
        }

        double slope_limited[5];
        for(int i=0; i<5; i++)
        {
            if ((W_plus[i]-W_equal[i])*(W_equal[i]-W_minus[i]) > 0 ){

                double slope_sign = 0.5*(W_plus[i]-W_minus[i])/ abs(0.5*(W_plus[i]-W_minus[i]));
                slope_limited[i] = min( abs(0.5*(W_plus[i]-W_minus[i])), 2*abs(W_equal[i]-W_minus[i]));
                slope_limited[i] = min( slope_limited[i], 2*abs(W_plus[i]-W_equal[i]))*slope_sign;
            }
            else{
                slope_limited[i] = 0;
            }
        }

        // get the face-centered variables
        for(int i=0; i<5; i++)
        {
            LR.L[pos][i] = W[pos][i] + 1/6*(W[pos][i]-W[pos-1][i])+1/3*(W[pos+1][i]-W[pos][i]);
            LR.R[pos][i] = W[pos][i] - 1/6*(W[pos+2][i]-W[pos+1][i])-1/3*(W[pos+1][i]-W[pos][i]);
            
            //ensure face-centered variable lie between nearby volume-average
            //discontinuity jump
            //double del2a = (W[pos+1][i]-2*W[pos][i]+W[pos-1][i])/(6*dx*dx);


        }

        // primitve variables -> conserved variables
        double L_ttemp[5] = {0};
        double R_ttemp[5] = {1};

        for(int i=0; i<5; i++) 
        {
            L_ttemp[i] = LR.L[pos][i];
            R_ttemp[i] = LR.R[pos][i];
        }

        vec5 L_temp = Primitive2Conserved(L_ttemp);
        vec5 R_temp = Primitive2Conserved(R_ttemp);
        
        for(int i=0; i<5; i++)
        {
            LR.L[pos][i] = L_temp.u[i];
            LR.R[pos][i] = R_temp.u[i];
            //printf("LR.L=%f\n",LR.L[pos][i]);
            //printf("LR.R=%f\n",LR.R[pos][i]); //Promblems before here.
        }
    }
    return LR;
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

int main()
{
    FILE *output;
    output = fopen("output.txt","w");

    omp_set_num_threads(N_thread);
    time_start = clock();
    double x;
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
    
    while (t >= 0)
    {
        // Set boundary condition
        BoundaryCondition( U );

        // calculate time-step
        dt = ComputeTimestep (U);
        //printf("dt = %f",dt);
        fprintf( output, "t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);
        printf("t = %13.7e --> %13.7e, dt = %13.7e\n", t, t+dt, dt);

        vecLR LR;
        //LR = DataReconstruction_PLM( U );
        LR = DataReconstruction_PPM( U );

        //updat the face-center variales by 0.5*dt
        for(int pos=1; pos<N-1; pos++)
        {
            double L_temp[5], R_temp[5];
            for(int i=0; i<5; i++)
            {
                L_temp[i] = LR.L[pos][i];
                R_temp[i] = LR.R[pos][i];
                //printf("LR.L=%f\n",LR.L[pos][i]);
            }

            vec5 Flux_L, Flux_R;
            Flux_L = Conserved2Flux(L_temp);
            Flux_R = Conserved2Flux(R_temp);
            
            double dflux[5] = {0};
            for(int i=0; i<5; i++)
            {
            dflux[i] = 0.5*dt/dx*(Flux_R.u[i]-Flux_L.u[i]);
            //printf("dflux[%d][%d]=%f\n",pos,i,dflux[i]);
            LR.L[pos][i] -= dflux[i];
            LR.R[pos][i] -= dflux[i];
            }
        }

        //printf("debug : compute fluxes\n");
        //compute fluxes
        double flux[N][5];
        # pragma omp parallel
        {
            # pragma omp for
            for(int pos=nghost; pos<N-nghost+1; pos++)
            {
                double R_temp[5], L_temp[5];
                for(int i=0; i<5; i++)
                {
                    R_temp[i] = LR.R[pos-1][i];
                    L_temp[i] = LR.L[pos][i];
                }
                // R_temp is the LEFT state at the j+1/2 interface
                vec5 flux_temp = {0};
                flux_temp = Roe(R_temp, L_temp);
                for(int i=0; i<5; i++)
                {
                    flux[pos][i] = flux_temp.u[i];
                    //printf("flux[%d][%d]=%f\n",pos,i,flux[pos][i]);
                }
            }
        }
        // update the volume-averaged input variables by dt
        # pragma omp parallel
        {
            # pragma omp for
            for(int pos= nghost; pos<N-nghost; pos++)
            {
                for(int i=0; i<5; i++)
                {
                    U[pos][i] -= dt/dx*(flux[pos+1][i]-flux[pos][i]);
                }
            }
        }
        
        // THIS LOOP IS FOR PRINT TO DATA
        for(int pos = nghost; pos<N-nghost; pos++)
        {
            double data_conservative[5] = {U[pos][0],U[pos][1],U[pos][2],U[pos][3],U[pos][4]};
            vec5 data_primitive = Conserved2Primitive(data_conservative);
        
            for(int i=0; i<5; i++)
                fprintf(output, "%.9e ",data_primitive.u[i]);
            fprintf(output, "\n");

        }        
        // update time
        t = t + dt;
        if(t >= end_time)
            break;
    }
    time_end = clock();
    printf("time = %d (ms)\n", time_end-time_start);
    fclose(output);
}
