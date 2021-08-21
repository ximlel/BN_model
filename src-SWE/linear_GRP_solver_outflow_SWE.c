#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"


double outflow_Riemann_solver_exact(double * u_star, double h_star, const double g, double u_L, double h_L, int * CRW, double tol, int N)
{
    double c_L = sqrt(g*h_L);
    double c_star = sqrt(g*h_star);
    double delta_h;

    //=====find out the kinds of the 1-wave
    if(h_star > h_L) // (u_star,h_star) lies on the shock branch of I1
	{
	    delta_h = h_star - h_L;
	    *u_star = delta_h * sqrt(0.5*g*(1.0/h_star + 1.0/h_L));
	    *u_star = u_L - (*u_star);
        CRW[0] = 0;
	}
    else // (u_star,h_star) lies on the rarefaction branch of I1
	{
	    *u_star = 2.0 * (c_star - c_L);
	    *u_star = u_L - (*u_star);
        CRW[0] = 1;
	}

    return 0.0;
}


void linear_GRP_solver_outflow_SWE(double * dire, double * mid, double * input)
{
    double h_star, h_L, s_h_L, u_L, s_u_L, s_u_R, B, s_B_L, s_B_R, g, eps;
    h_star = input[0];
    h_L = input[1];
    s_h_L = input[2];
    u_L = input[3];
    s_u_L = input[4];
    B = input[5];
    s_B_L = input[6];
    s_B_R = input[7];
    g = input[8];
    eps = input[9];

    double c_L, c_R;
    c_L = sqrt(g * h_L);

    int CRW[2];
    double u_star, c_star;
    double a_L, b_L, d_L;
    double pt_pa, ps_pa, pt_pb, pr_pb, alpha, beta, phi, det_A;
    double s_L, s_s_L;
    s_L = u_L + 2.0*c_L;
    s_s_L = s_u_L + sqrt(g/h_L)*s_h_L;
    double d_G_L[2], d_G[2];
    double A_gI_L[2][2], A_gI[2][2];
    double speed_L;
    double ss_shk_spd, dh_dt, du_dt;

    //=========non-acoustic case==========
    outflow_Riemann_solver_exact(&u_star, h_star, g, u_L, h_L, CRW, eps, 50);
    c_star = sqrt(g*h_star);
    
    //----------solving the LINEAR GRP----------
    if(CRW[0])
	speed_L = u_L - c_L;
    else
	speed_L = u_L - h_star * sqrt(0.5*g * (1.0/h_star + 1.0/h_L));

    //------trivial case------
    if(speed_L > eps) //the t-axe is on the left side of all the three waves
	{
	    dire[0] = -s_h_L*u_L - h_L*s_u_L;
	    dire[1] = -s_h_L - u_L*s_u_L/g - s_B_L;
	    dire[1]*= g;

	    mid[0] = h_L;
	    mid[1] = u_L;
	}
    //----non-trivial case----
    else
	{
	    if((CRW[0]) && ((u_star-c_star) > 0.0)) // the t-axe is in a 1-CRW
		{
		    mid[1] = (u_L+2.0*c_L)/3.0;
		    mid[0] = mid[1]*mid[1]/g;
		    phi = -g*s_B_L;
		    a_L = sqrt(g/mid[0]);
		    b_L = 1.0;
		    beta = 0.0;
		    ps_pa = (phi-2.0*c_L*s_s_L)/pow(3.0*c_L,1.5);
		    ps_pa-= phi*(pow(s_L-beta,-1.5)-pow(3.0*c_L,-1.5));
		    pt_pa = 1.0/pow(s_L-beta,1.5);
		    d_L = ps_pa/pt_pa;
		}
	    //--non-sonic case--
	    else
		{
		    //determine a_L, b_L and d_L
		    if(CRW[0]) //the 1-wave is a CRW
			{
			    phi = -g*s_B_L;
			    det_A = u_star*u_star - g*h_star;
			    a_L = ( sqrt(g/h_star)*u_star - g)/det_A;
			    b_L = (-sqrt(g/h_star)*h_star + u_star)/det_A;
			    beta = u_star - c_star;
			    ps_pa = (phi-2.0*c_L*s_s_L)/pow(3.0*c_L,1.5);
			    ps_pa-= phi*(pow(3.0*c_star,-1.5)-pow(3.0*c_L,-1.5));
			    pt_pa = 1.0/pow(s_L-beta,1.5);
			    d_L = ps_pa/pt_pa - (u_star-c_star)*phi*b_L;
			    a_L*= 2.0*c_star;
			    b_L*= 2.0*c_star;
			}
		    else //the 1-wave is a shock
			{
			    phi = -g*s_B_L;
			    det_A = u_star*u_star - g*h_star;
			    d_G_L[0] = (h_L-h_star)*(-0.5*g/h_L   *(2.0*h_L*h_L+h_star*h_star+h_L*h_star));
			    d_G_L[1] = (h_L-h_star)*(2.0*h_L*h_star*(u_L-u_star)/(h_L-h_star));
			    d_G[0]   = (h_star-h_L)*(-0.5*g/h_star*(2.0*h_star*h_star+h_L*h_L+h_L*h_star));
			    d_G[1]   = (h_star-h_L)*(2.0*h_L*h_star*(u_L-u_star)/(h_L-h_star));
			    A_gI_L[0][0] = sqrt(0.5*g*(h_L+h_star)*h_star/h_L);
			    A_gI_L[0][1] = h_L;
			    A_gI_L[1][0] = g;
			    A_gI_L[1][1] = sqrt(0.5*g*(h_L+h_star)*h_star/h_L);
			    A_gI[0][0]   = sqrt(0.5*g*(h_L+h_star)*h_L/h_star);
			    A_gI[0][1]   = h_star;
			    A_gI[1][0]   = g;
			    A_gI[1][1]   = sqrt(0.5*g*(h_L+h_star)*h_L/h_star);
			    a_L = (d_G[0]*A_gI[0][0]+d_G[1]*A_gI[1][0])*u_star - (d_G[0]*A_gI[0][1]+d_G[1]*A_gI[1][1])*g;
			    a_L/= det_A;
			    b_L =-(d_G[0]*A_gI[0][0]+d_G[1]*A_gI[1][0])*h_star + (d_G[0]*A_gI[0][1]+d_G[1]*A_gI[1][1])*u_star;
			    b_L/= det_A;
			    d_L = (d_G_L[0]*A_gI_L[0][0]+d_G_L[1]*A_gI_L[1][0])*s_h_L + (d_G_L[0]*A_gI_L[0][1]+d_G_L[1]*A_gI_L[1][1])*s_u_L;
			    d_L-= speed_L*phi*(-h_star*d_G[0]+u_star*d_G[1])/det_A;
			    d_L-= phi*d_G_L[1];
			}
		    mid[0] = h_star;
		    mid[1] = u_star;
		    dire[0] = 0.0;
		    dire[1] = d_L/b_L;
		    //--end of non-sonic case--
		}
	    // stationary shocks case
	    if ((!CRW[0]) && fabs(speed_L) <= eps)
		{
		    phi = -g*s_B_L;
		    dh_dt = -s_h_L*u_L - h_L*s_u_L;
		    du_dt = -s_h_L - u_L*s_u_L/g - s_B_L;
		    du_dt*= g;
		    ss_shk_spd = 1.0/(h_star - h_L)*(u_star*dire[0] + h_star*dire[1] - u_L*dh_dt - h_L*du_dt);
		    if (ss_shk_spd > 0.0)
			{
			    mid[0] = h_L;
			    mid[1] = u_L;
			    dire[0] = dh_dt;
			    dire[1] = du_dt;		    
			}
		}
	    //----end of non-trivial case----
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    double *M;
    int m,n;
    M = mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
     
    plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);
    double *A;
    A = mxGetPr(plhs[0]);
     
    linear_GRP_solver_outflow_SWE(A,A+2,M);
}
