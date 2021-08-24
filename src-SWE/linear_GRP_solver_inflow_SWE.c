#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"


double inflow_Riemann_solver_exact(double * u_star, double * h_star, const double g, double q_star, double u_R, double h_R, int * CRW, double tol, int N)
{
    double c_R = sqrt(g*h_R);
    double delta_h;
    
    double k1, k3, h_INT, h_INT0, u_INT;
    double v_L, v_R, gap;
    double temp1, temp2;
    int n = 0;

    //=====find out the kind of the 3-wave
    if(q_star > h_R*u_R) // (u_star, h_star) lies on the shock branch of I3       
	CRW[1] = 0;
    else // (u_star, h_star) lies on the rarefaction branch of I3
	CRW[1] = 1;
	
    *h_star = h_R;
    *u_star = q_star/h_R;        

    //======one step of the Newton ietration to get the intersection point of   and I3====
    k1 = -q_star/(*h_star)/(*h_star);//the (h,u)-tangent slope on I1 at (u_star,h_star), i.e. [du/dh](h_star)

    temp1 = 1.0/(*h_star) + 2.0/h_R + h_R/(*h_star)/(*h_star);
    temp2 = 2.0*sqrt(1.0/(*h_star)+1.0/h_R);
    k3 =  sqrt(0.5 * g) * temp1 / temp2;//the (h,u)-tangent slope on I3 at (u_R,h_R), i.e. [du/dh](h_R)
    //the intersect of (u-u_star)=k1*(h-h_star) and (u-u_R)=k3*(h-h_R)
    h_INT = h_R + (u_R - (*u_star)) / (k1 - k3);    
    if(h_INT <= 0)
	h_INT = h_R;
    h_INT = 0.5*h_INT;
    
    //=======compute the gap between U^n_star and U^n_L(see Appendix C)=======
    v_L = q_star/h_INT;
    
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
	    k1 = -q_star/h_INT/h_INT;
		
	    //the (h,u)-tangent slope on I3 at (v_R,h_INT), i.e. [du/dh](h_INT)	    
	    temp1 = 1.0/h_INT + 2.0/h_R + h_R/h_INT/h_INT;
	    temp2 = 2.0*sqrt(1.0/h_INT + 1.0/h_R);
	    k3 = sqrt(0.5 * g) * temp1 / temp2;

	    //the intersect of (u-v_L)=k1*(h-h_INT) and (u-v_R)=k3*(h-h_INT)
	    h_INT0 = h_INT + (v_R - v_L) / (k1 - k3);
	    if(h_INT0 < 0.0)
		h_INT = 0.5*h_INT;
	    else
		h_INT = h_INT0;

	    //------the gap------
	    ++n;
	    v_L = q_star/h_INT;
        
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


void linear_GRP_solver_inflow_SWE(double * dire, double * mid, double * input)
{
    double q_star, h_R, s_h_R, u_R, s_u_R, B, s_B_L, s_B_R, g, eps;
    q_star = input[0];
    h_R = input[1];
    s_h_R = input[2];
    u_R = input[3];
    s_u_R = input[4];
    B = input[5];
    s_B_L = input[6];
    s_B_R = input[7];
    g = input[8];
    eps = input[9];

    double c_R;
    c_R = sqrt(g * h_R);

    int CRW[2];
    double u_star, h_star, c_star;
    double a_R, b_R, d_R;
    double pt_pa, ps_pa, pt_pb, pr_pb, alpha, beta, phi, det_A;
    double r_R, s_r_R;
    r_R = u_R - 2.0*c_R;
    s_r_R = s_u_R - sqrt(g/h_R)*s_h_R;
    double d_G_R[2], d_G[2];
    double A_gI_R[2][2], A_gI[2][2];
    double speed_R;
    double ss_shk_spd, dh_dt, du_dt;

    //=========non-acoustic case==========
    inflow_Riemann_solver_exact(&u_star, &h_star, g, q_star, u_R, h_R, CRW, eps, 50);
    c_star = sqrt(g*h_star);
    
    //----------solving the LINEAR GRP----------
    if(CRW[1])
	speed_R = u_R + c_R;
    else
	speed_R = u_R + h_star * sqrt(0.5*g * (1.0/h_star + 1.0/h_R));

    //------trivial case------
    if(speed_R < -eps) //the t-axe is on the right side of all the three waves
	{
	    dire[0] = -s_h_R*u_R - h_R*s_u_R;
	    dire[1] = -s_h_R - u_R*s_u_R/g - s_B_R;
	    dire[1]*= g;

	    mid[0] = h_R;
	    mid[1] = u_R;
	}
    else
	{
	    if((CRW[1]) && ((u_star+c_star) < 0.0)) // the t-axe is in a 3-CRW
		{
		    mid[1] = (u_R-2.0*c_R)/3.0;
		    mid[0] = mid[1]*mid[1]/g;
		    phi = -g*s_B_R;
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
		    dire[0] = (- d_R*u_star)/(u_star*b_R - a_R*h_star);
		    dire[1] = (- d_R*h_star)/(a_R*h_star - u_star*b_R);
		    //--end of non-sonic case--
		}
	    // stationary shocks case
	    if ((!CRW[1]) && fabs(speed_R) <= eps)
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
     
    linear_GRP_solver_inflow_SWE(A,A+2,M);
}
