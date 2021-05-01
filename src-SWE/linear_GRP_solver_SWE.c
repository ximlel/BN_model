#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"


double Riemann_solver_exact(double * u_star, double * h_star, const double g, double u_L, double u_R, double h_L, double h_R, int * CRW, double tol, int N)
{
    double c_L = sqrt(g*h_L);
    double c_R = sqrt(g*h_R);
    double delta_h, u_LR, u_RL;
    
    double k1, k3, h_INT, h_INT0, u_INT;
    double v_L, v_R, gap;
    double temp1, temp2;
    int n = 0;

    //=====find out the kinds of the 1-wave and the 3-wave
    if(h_R > h_L) // (u_LR,h_R) lies on the shock branch of I1
	{
	    delta_h = h_R - h_L;
	    u_LR = delta_h * sqrt(0.5*g*(1.0/h_R + 1.0/h_L));
	    u_LR = u_L - u_LR;
	}
    else // (u_LR,h_R) lies on the rarefaction branch of I1
	{
	    u_LR = 2.0 * (c_R - c_L);
	    u_LR = u_L - u_LR;
	}
    if(h_L > h_R) // (u_RL, h_L) lies on the shock branch of I3
	{
	    delta_h = h_L - h_R;
	    u_RL = delta_h * sqrt(0.5*g*(1.0/h_L + 1.0/h_R));
	    u_RL = u_R + u_RL;
	}
    else // (u_RL, h_L) lies on the rarefaction branch of I3
	{
	    u_RL = 2.0 * (c_L - c_R);
	    u_RL = u_R + u_RL;
	}
    if(u_LR > u_R+tol)
	CRW[1] = 0;
    else
	CRW[1] = 1;
    if(u_RL > u_L-tol)
	CRW[0] = 1;
    else
	CRW[0] = 0;

    //======one step of the Newton ietration to get the intersection point of I1 and I3====
    k1 = -sqrt(g/h_L);//the (h,u)-tangent slope on I1 at (u_L,h_L), i.e. [du/dh](h_L)
    k3 =  sqrt(g/h_R);//the (h,u)-tangent slope on I3 at (u_R,h_R), i.e. [du/dh](h_R)
    //the intersect of (u-u_L)=k1*(h-h_L) and (u-u_R)=k3*(h-h_R)
    h_INT = (k1*h_L - k3*h_R - u_L + u_R) / (k1 - k3);
    if(h_INT <= 0)
	h_INT = (h_L<h_R)? h_L : h_R;
    h_INT = 0.5*h_INT;
    
    //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
    if(h_INT > h_L)
	{
	    delta_h = h_INT - h_L;
	    v_L = delta_h * sqrt(0.5*g*(1.0/h_INT + 1.0/h_L));
	    v_L = u_L - v_L;
	}
    else
	{
	    v_L = 2.0 * (sqrt(g*h_INT) - c_L);
	    v_L = u_L - v_L;
	}
    if(h_INT > h_R)
	{
	    delta_h = h_INT - h_R;
	    v_R = delta_h * sqrt(0.5*g*(1.0/h_INT + 1.0/h_R));
	    v_R = u_R + v_R;
	}
    else
	{    
	    v_R = 2.0 * (sqrt(g*h_INT) - c_R);
	    v_R = u_R + v_R;
	}
    gap = fabs(v_L - v_R);


    //=======THE NEWTON ITERATION=======
    while((gap > tol) && (n != N))
	{
	    //the (h,u)-tangent slope on I1 at (v_L,h_INT), i.e. [du/dh](h_INT)
	    if(h_INT > h_L)
		{
		    temp1 = 1.0/h_INT + 2.0/h_L + h_L/h_INT/h_INT;
		    temp2 = 2.0*sqrt(1.0/h_INT + 1.0/h_L);
		    k1 = -sqrt(0.5 * g) * temp1 / temp2;
		}
	    else
		{
		    k1 = -sqrt(g/h_INT);
		}
	    //the (h,u)-tangent slope on I3 at (v_R,h_INT), i.e. [du/dh](h_INT)
	    if(h_INT > h_R)
		{
		    temp1 = 1.0/h_INT + 2.0/h_R + h_R/h_INT/h_INT;
		    temp2 = 2.0*sqrt(1.0/h_INT + 1.0/h_R);
		    k3 = sqrt(0.5 * g) * temp1 / temp2;
		}
	    else
		{
		    k3 = sqrt(g/h_INT);
		}

	    //the intersect of (u-v_L)=k1*(h-h_INT) and (u-v_R)=k3*(h-h_INT)
	    h_INT0 = h_INT + (v_R - v_L) / (k1 - k3);
	    if(h_INT0 < 0.0)
		h_INT = 0.5*h_INT;
	    else
		h_INT = h_INT0;

	    //------the gap------
	    ++n;
	    if(h_INT > h_L)
		{
		    delta_h = h_INT - h_L;
		    v_L = delta_h * sqrt(0.5*g*(1.0/h_INT + 1.0/h_L));
		    v_L = u_L - v_L;
		}
	    else
		{
		    v_L = 2.0 * (sqrt(g*h_INT) - c_L);
		    v_L = u_L - v_L;
		}
	    if(h_INT > h_R)
		{
		    delta_h = h_INT - h_R;
		    v_R = delta_h * sqrt(0.5*g*(1.0/h_INT + 1.0/h_R));
		    v_R = u_R + v_R;
		}
	    else
		{    
		    v_R = 2.0 * (sqrt(g*h_INT) - c_R);
		    v_R = u_R + v_R;
		}
	    gap = fabs(v_L - v_R);
	}
    u_INT = k1*(v_R-v_L)/(k1-k3) + v_L;

    *h_star = h_INT;
    *u_star = u_INT;

    return gap;
}


void linear_GRP_solver_SWE(double * dire, double * mid, double * input)
{
    double h_L, h_R, s_h_L, s_h_R, u_L, u_R, s_u_L, s_u_R, B, s_B_L, s_B_R, g, eps;
    h_L = input[0];
    h_R = input[1];
    s_h_L = input[2];
    s_h_R = input[3];
    u_L = input[4];
    u_R = input[5];
    s_u_L = input[6];
    s_u_R = input[7];
    B = input[8];
    s_B_L = input[9];
    s_B_R = input[10];
    g = input[11];
    eps = input[12];

    double c_L, c_R;
    c_L = sqrt(g * h_L);
    c_R = sqrt(g * h_R);

    int CRW[2];
    double u_star, h_star, c_star;
    double a_L, b_L, d_L, a_R, b_R, d_R;
    double pt_pa, ps_pa, pt_pb, pr_pb, alpha, beta, phi, det_A;
    double s_L, s_s_L, r_R, s_r_R;
    s_L = u_L + 2.0*c_L;
    s_s_L = s_u_L + sqrt(g/h_L)*s_h_L;
    r_R = u_R - 2.0*c_R;
    s_r_R = s_u_R - sqrt(g/h_R)*s_h_R;
    double d_G_L[2], d_G_R[2], d_G[2];
    double A_gI_L[2][2], A_gI_R[2][2], A_gI[2][2];
    double speed_L, speed_R;
    double ss_shk_spd, dh_dt, du_dt;

    //=========non-acoustic case==========
    Riemann_solver_exact(&u_star, &h_star, g, u_L, u_R, h_L, h_R, CRW, eps, 50);
    c_star = sqrt(g*h_star);
    
    //----------solving the LINEAR GRP----------
    if(CRW[0])
	speed_L = u_L - c_L;
    else
	speed_L = u_L - h_star * sqrt(0.5*g * (1.0/h_star + 1.0/h_L));
    if(CRW[1])
	speed_R = u_R + c_R;
    else
	speed_R = u_R + h_star * sqrt(0.5*g * (1.0/h_star + 1.0/h_R));

    //------trivial case------
    if(speed_L > eps) //the t-axe is on the left side of all the three waves
	{
	    dire[0] = -s_h_L*u_L - h_L*s_u_L;
	    dire[1] = -s_h_L - u_L*s_u_L/g - s_B_L;
	    dire[1]*= g;

	    mid[0] = h_L;
	    mid[1] = u_L;
	}
    else if(speed_R < -eps) //the t-axe is on the right side of all the three waves
	{
	    dire[0] = -s_h_R*u_R - h_R*s_u_R;
	    dire[1] = -s_h_R - u_R*s_u_R/g - s_B_R;
	    dire[1]*= g;

	    mid[0] = h_R;
	    mid[1] = u_R;
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
		    a_R = -sqrt(g/mid[0]);
		    b_R = 1.0;
		    d_R = phi;
		}
	    else if((CRW[1]) && ((u_star+c_star) < 0.0)) // the t-axe is in a 3-CRW
		{
		    mid[1] = (u_R-2.0*c_R)/3.0;
		    mid[0] = mid[1]*mid[1]/g;
		    phi = -g*s_B_R;
		    a_L = sqrt(g/mid[0]);
		    b_L = 1.0;
		    d_L = phi;
		    a_R = -sqrt(g/mid[0]);
		    b_R = 1.0;
		    alpha = 0.0;
		    pr_pb = (phi+2.0*c_R*s_r_R)/pow(3*c_R,1.5);
		    pr_pb+= phi*(pow(alpha-r_R,-1.5)-pow(3.0*c_R,-1.5));
		    pt_pb = 1.0/pow(alpha-r_R,1.5);
		    d_R = pr_pb/pt_pb;		    
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
		    //determine a_R, b_R and d_R
		    if(CRW[1]) //the 3-wave is a CRW
			{
			    phi = -g*s_B_R;
			    det_A = u_star*u_star - g*h_star;
			    a_R = (-sqrt(g/h_star)*u_star - g)/det_A;
			    b_R = ( sqrt(g/h_star)*h_star + u_star)/det_A;
			    alpha = u_star + c_star;
			    pr_pb = (phi+2.0*c_R*s_r_R)/pow(3*c_R,1.5);
			    pr_pb = pr_pb + phi*(pow(3.0*c_star,-1.5)-pow(3.0*c_R,-1.5));
			    pt_pb = 1.0/pow(alpha-r_R,1.5);
			    d_R = pr_pb/pt_pb - (u_star+c_star)*phi*b_R;	
			    a_R*=-2.0*c_star;
			    b_R*=-2.0*c_star;
			}
		    else //the 3-wave is a shock
			{
			    phi = -g*s_B_R;
			    det_A = u_star*u_star - g*h_star;
			    d_G_R[0] = (h_R-h_star)*(-0.5*g/h_R   *(2.0*h_R*h_R+h_star*h_star+h_R*h_star));
			    d_G_R[1] = (h_R-h_star)*(2.0*h_R*h_star*(u_R-u_star)/(h_R-h_star));
			    d_G[0]   = (h_star-h_R)*(-0.5*g/h_star*(2.0*h_star*h_star+h_R*h_R+h_R*h_star));
			    d_G[1]   = (h_star-h_R)*(2.0*h_R*h_star*(u_R-u_star)/(h_R-h_star));
			    A_gI_R[0][0] = -sqrt(0.5*g*(h_R+h_star)*h_star/h_R);
			    A_gI_R[0][1] = h_R;
			    A_gI_R[1][0] = g;
			    A_gI_R[1][1] = -sqrt(0.5*g*(h_R+h_star)*h_star/h_R);
			    A_gI[0][0]   = -sqrt(0.5*g*(h_R+h_star)*h_R/h_star);
			    A_gI[0][1]   = h_star;
			    A_gI[1][0]   = g;
			    A_gI[1][1]   = -sqrt(0.5*g*(h_R+h_star)*h_R/h_star);
			    a_R = (d_G[0]*A_gI[0][0]+d_G[1]*A_gI[1][0])*u_star - (d_G[0]*A_gI[0][1]+d_G[1]*A_gI[1][1])*g;
			    a_R/= det_A;
			    b_R =-(d_G[0]*A_gI[0][0]+d_G[1]*A_gI[1][0])*h_star + (d_G[0]*A_gI[0][1]+d_G[1]*A_gI[1][1])*u_star;
			    b_R/= det_A;
			    d_R = (d_G_R[0]*A_gI_R[0][0]+d_G_R[1]*A_gI_R[1][0])*s_h_R + (d_G_R[0]*A_gI_R[0][1]+d_G_R[1]*A_gI_R[1][1])*s_u_R;
			    d_R-= speed_R*phi*(-h_star*d_G[0]+u_star*d_G[1])/det_A;
			    d_R-= phi*d_G_R[1];
			}
		    mid[0] = h_star;
		    mid[1] = u_star;
		    dire[0] = (d_L*b_R - d_R*b_L)/(a_L*b_R - a_R*b_L);
		    dire[1] = (d_L*a_R - d_R*a_L)/(a_R*b_L - a_L*b_R);
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
	    else if ((!CRW[1]) && fabs(speed_R) <= eps)
		{
		    phi = -g*s_B_R;
		    dh_dt = -s_h_R*u_R - h_R*s_u_R;
		    du_dt = -s_h_R - u_R*s_u_R/g - s_B_R;
		    du_dt*= g;
		    ss_shk_spd = 1.0/(h_star - h_R)*(u_star*dire[0] + h_star*dire[1] - u_R*dh_dt - h_R*du_dt);
		    if (ss_shk_spd > 0.0)
			{
			    mid[0] = h_R;
			    mid[1] = u_R;
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
     
    linear_GRP_solver_SWE(A,A+2,M);
}
