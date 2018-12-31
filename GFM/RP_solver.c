#include <math.h>
#include <stdio.h>
#include "mex.h"


double RP_solver(double * star, double *input)
{
  double gammaL, gammaR, u_L, u_R, p_L, p_R, rho_L, rho_R, c_L, c_R, eps, tol;
  int N=500;
  rho_L=input[0];
  rho_R=input[1];
  u_L=input[2];
  u_R=input[3];
  p_L=input[4];
  p_R=input[5];
  gammaL=input[6];
  gammaR=input[7];
  eps=input[8];
  tol=eps;
  c_L=sqrt(gammaL*p_L/rho_L);
  c_R=sqrt(gammaR*p_R/rho_R);

  //double zeta_l, zeta_r;
  double muL, nuL, sigmaL, gammaL_1;
  double muR, nuR, sigmaR, gammaR_1;
  //double c_L, c_R;
  double delta_p, u_LR, u_RL;
  double k1, k3, p_INT, p_INT0, u_INT;
  double rho_star_l, rho_star_r;
  double v_L, v_R, gap;
  double temp1, temp2, temp3;
  double dbg = 2;
  int n = 0;

  muL = (gammaL-1.0) / (2.0*gammaL);
  nuL = (gammaL+1.0) / (2.0*gammaL);
  sigmaL = (gammaL - 1.0) / (gammaL + 1.0);
  gammaL_1 = gammaL-1.0;
  muR = (gammaR-1.0) / (2.0*gammaR);
  nuR = (gammaR+1.0) / (2.0*gammaR);
  sigmaR = (gammaR - 1.0) / (gammaR + 1.0);
  gammaR_1 = gammaR-1.0;

  //c_L = sqrt(gamma * p_L / rho_l);
  //c_R = sqrt(gamma * p_R / rho_r);
  //rho_l = 1.0 / rho_l;
  //rho_r = 1.0 / rho_r;

  //zeta_l = pow(p_L, (gamma-1.0) / (2.0*gamma));
  //zeta_r = pow(p_R, (gamma-1.0) / (2.0*gamma));

  //=====find out the kinds of the 1-wave and the 3-wave, page 132 in the GRP book
  //find out where (u_LR,p_R) lies on the curve of LEFT state
  // (u_LR,p_R) lies on the shock branch of I1
  {
    delta_p = p_R - p_L;
    u_LR = sqrt(1.0 + nuL*delta_p/p_L);
    u_LR = delta_p * c_L / gammaL / p_L / u_LR;
    u_LR = u_L - u_LR;
  }
  //find out where (u_RL,p_L) lies on the curve of RIGHT state
  // (u_RL, p_L) lies on the shock branch of I3
  {
    delta_p = p_L - p_R;
    u_RL = sqrt(1.0 + nuR*delta_p/p_R);
    u_RL = delta_p * c_R / gammaR / p_R / u_RL;
    u_RL = u_R + u_RL;
  }
  
  //======one step of the Newton ietration to get the intersection point of I1 and I3====
  k1 = -c_L / p_L / gammaL;//the (p,u)-tangent slope on I1 at (u_L,p_L), i.e. [du/dp](p_L)
  k3 =  c_R / p_R / gammaR;//the (p,u)-tangent slope on I3 at (u_R,p_R), i.e. [du/dp](p_R)
  //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
  p_INT = (k1*p_L - k3*p_R - u_L + u_R) / (k1 - k3);
  if(p_INT < 0)
    p_INT = (p_L<p_R)? p_L : p_R;
  p_INT = 0.5*p_INT;

  //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
  {
    delta_p = p_INT - p_L;
    v_L = sqrt(1.0 + nuL*delta_p/p_L);
    v_L = delta_p * c_L / gammaL / p_L / v_L;
    v_L = u_L - v_L;
  }
  {
    delta_p = p_INT - p_R;
    v_R = sqrt(1.0 + nuR*delta_p/p_R);
    v_R = delta_p * c_R / gammaR / p_R / v_R;
    v_R = u_R + v_R;
  }

  gap = fabs(v_L - v_R);


  //=======THE NEWTON ITERATION=======
  while((gap > tol) && (n != N))
  {
    //the (p,u)-tangent slope on I1 at (v_L,p_INT), i.e. [du/dp](p_INT)
    {
      delta_p = p_INT - p_L;
      temp1 = 1.0 / sqrt(1.0 + nuL*delta_p/p_L);
      temp2 = c_L / gammaL / p_L;
      temp3 = 0.5 * temp2 * nuL / p_L;
      k1 = temp3*delta_p*pow(temp1,3.0) - temp2*temp1;
    }
    //the (p,u)-tangent slope on I3 at (v_R,p_INT), i.e. [du/dp](p_INT)
    {
      delta_p = p_INT - p_R;
      temp1 = 1.0 / sqrt(1.0 + nuR*delta_p/p_R);
      temp2 = c_R / gammaR / p_R;
      temp3 = 0.5 * temp2 * nuR / p_R;
      k3 = temp2*temp1 - temp3*delta_p*pow(temp1,3.0);
    }

    //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
    p_INT0 = p_INT + (v_R - v_L) / (k1 - k3);
    if(p_INT0 < 0.0)
      p_INT = 0.5*p_INT;
    else
      p_INT = p_INT0;

    //------the gap------
    ++n;
    {
      delta_p = p_INT - p_L;
      v_L = sqrt(1.0 + nuL*delta_p/p_L);
      v_L = delta_p * c_L / gammaL / p_L / v_L;
      v_L = u_L - v_L;
    }
    {
      delta_p = p_INT - p_R;
      v_R = sqrt(1.0 + nuR*delta_p/p_R);
      v_R = delta_p * c_R / gammaR / p_R / v_R;
      v_R = u_R + v_R;
    }

    gap = fabs(v_L - v_R);
  }
//printf("n=%d\tN=%d\n", n, N);

  u_INT = k1*(v_R-v_L)/(k1-k3)+v_L;

  star[0] = p_INT;
  star[1] = u_INT;
  const double zetaL = (gammaL-1.0)/(gammaL+1.0);
  const double zetaR = (gammaR-1.0)/(gammaR+1.0);
  star[2] = rho_L*(p_INT+zetaL*p_L)/(p_L+zetaL*p_INT);
  star[3] = rho_R*(p_INT+zetaR*p_R)/(p_R+zetaR*p_INT);

  return gap;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    double *M;
    int m,n;
    M = mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
     
    plhs[0] = mxCreateDoubleMatrix(9,1,mxREAL);
    double *A;
    A = mxGetPr(plhs[0]);
     
    RP_solver(A,M);
}
