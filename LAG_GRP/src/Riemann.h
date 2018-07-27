#ifndef _RIEMANN_H
#define _RIEMANN_H

double Riemann_solver_exact(double &P_star, double &U_star, double &rho_starL,double &rho_starR,double rho_L,double rho_R, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, double gamma_L, double gamma_R)
{
  int CRW[2];
  double tol=EPS;
  const int N=500;

  double mu_L, nu_L, sigma_L, mu_R, nu_R, sigma_R;
  double delta_p, u_LR, u_RL;
  double k1, k3, p_INT, p_INT0, u_INT;
  double rho_star_l, rho_star_r;
  double v_L, v_R, gap;
  double temp1, temp2, temp3;
  double dbg = 2;
  int n = 0;

  mu_L = (gamma_L-1.0) / (2.0*gamma_L);
  nu_L = (gamma_L+1.0) / (2.0*gamma_L);
  sigma_L = (gamma_L - 1.0) / (gamma_L + 1.0);
  mu_R = (gamma_R-1.0) / (2.0*gamma_R);
  nu_R = (gamma_R+1.0) / (2.0*gamma_R);
  sigma_R = (gamma_R - 1.0) / (gamma_R + 1.0);

  //=====find out the kinds of the 1-wave and the 3-wave
  if(p_R > p_L) // (u_LR,p_R) lies on the shock branch of I1
  {
    delta_p = p_R - p_L;
    u_LR = sqrt(1.0 + nu_L*delta_p/p_L);
    u_LR = delta_p * c_L / gamma_L / p_L / u_LR;
    u_LR = u_L - u_LR;
  }
  else // (u_LR,p_R) lies on the rarefaction branch of I1
  {
    u_LR = pow(p_R/p_L, mu_L) - 1.0;
    u_LR = 2.0 * c_L * u_LR / (gamma_L-1.0);
    u_LR = u_L - u_LR;
  }
  if(p_L > p_R) // (u_RL, p_L) lies on the shock branch of I3
  {
    delta_p = p_L - p_R;
    u_RL = sqrt(1.0 + nu_R*delta_p/p_R);
    u_RL = delta_p * c_R / gamma_R / p_R / u_RL;
    u_RL = u_R + u_RL;
  }
  else // (u_RL, p_L) lies on the rarefaction branch of I3
  {
    u_RL = pow(p_L/p_R, mu_R) - 1.0;
    u_RL = 2.0 * c_R * u_RL / (gamma_R-1.0);
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
  k1 = -c_L / p_L / gamma_L;//the (p,u)-tangent slope on I1 at (u_L,p_L), i.e. [du/dp](p_L)
  k3 =  c_R / p_R / gamma_R;//the (p,u)-tangent slope on I3 at (u_R,p_R), i.e. [du/dp](p_R)
  //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
  p_INT = (k1*p_L - k3*p_R - u_L + u_R) / (k1 - k3);
  if(p_INT < 0)
    p_INT = (p_L<p_R)? p_L : p_R;
  p_INT = 0.5*p_INT;

  //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
  if(p_INT > p_L)
  {
    delta_p = p_INT - p_L;
    v_L = sqrt(1.0 + nu_L*delta_p/p_L);
    v_L = delta_p * c_L / gamma_L / p_L / v_L;
    v_L = u_L - v_L;
  }
  else
  {
    v_L = pow(p_INT/p_L, mu_L) - 1.0;
    v_L = 2.0 * c_L * v_L / (gamma_L-1.0);
    v_L = u_L - v_L;
  }
  if(p_INT > p_R)
  {
    delta_p = p_INT - p_R;
    v_R = sqrt(1.0 + nu_R*delta_p/p_R);
    v_R = delta_p * c_R / gamma_R / p_R / v_R;
    v_R = u_R + v_R;
  }
  else
  {
    dbg = pow(p_INT/p_R, mu_R);
    v_R = pow(p_INT/p_R, mu_R) - 1.0;
    v_R = 2.0 * c_R * v_R / (gamma_R-1.0);
    v_R = u_R + v_R;
  }
  gap = fabs(v_L - v_R);
/*
if (fabs(u_L - u_R) < tol && fabs(p_L - p_R) < tol)
{
  P_star = 0.5*(p_L + p_R);
  U_star = 0.5*(u_L + u_R);

  return fabs(u_L - u_R);
}
*/
  //=======THE NEWTON ITERATION=======
  while((gap > tol) && (n != N))
  {
    //the (p,u)-tangent slope on I1 at (v_L,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_L)
    {
      delta_p = p_INT - p_L;
      temp1 = 1.0 / sqrt(1.0 + nu_L*delta_p/p_L);
      temp2 = c_L / gamma_L / p_L;
      temp3 = 0.5 * temp2 * nu_L / p_L;
      k1 = temp3*delta_p*pow(temp1,3.0) - temp2*temp1;
    }
    else
    {
      temp2 = c_L / gamma_L / p_L;
      temp1 = 1.0 / pow(p_INT/p_L, nu_L);
      k1 = -temp1 * temp2;
    }
    //the (p,u)-tangent slope on I3 at (v_R,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      temp1 = 1.0 / sqrt(1.0 + nu_R*delta_p/p_R);
      temp2 = c_R / gamma_R / p_R;
      temp3 = 0.5 * temp2 * nu_R / p_R;
      k3 = temp2*temp1 - temp3*delta_p*pow(temp1,3.0);
    }
    else
    {
      temp2 = c_R / gamma_R / p_R;
      temp1 = 1.0 / pow(p_INT/p_R, nu_R);
      k3 = temp1 * temp2;
    }

    //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
    p_INT0 = p_INT + (v_R - v_L) / (k1 - k3);
    if(p_INT0 < 0.0)
      p_INT = 0.5*p_INT;
    else
      p_INT = p_INT0;

    //------the gap------
    ++n;
    if(p_INT > p_L)
    {
      delta_p = p_INT - p_L;
      v_L = sqrt(1.0 + nu_L*delta_p/p_L);
      v_L = delta_p * c_L / gamma_L / p_L / v_L;
      v_L = u_L - v_L;
    }
    else
    {
      v_L = pow(p_INT/p_L, mu_L) - 1.0;
      v_L = 2.0 * c_L * v_L / (gamma_L-1.0);
      v_L = u_L - v_L;
    }
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      v_R = sqrt(1.0 + nu_R*delta_p/p_R);
      v_R = delta_p * c_R / gamma_R / p_R / v_R;
      v_R = u_R + v_R;
    }
    else
    {
      v_R = pow(p_INT/p_R, mu_R) - 1.0;
      v_R = 2.0 * c_R * v_R / (gamma_R-1.0);
      v_R = u_R + v_R;
    }

    gap = fabs(v_L - v_R);
  }


  u_INT = k1*(v_R-v_L)/(k1-k3)+v_L;

  P_star = p_INT;
  U_star = u_INT;

	if(P_star<=p_L)//Left rarefaction wave
		rho_starL=rho_L*pow(P_star/p_L,1./gamma_L);
	else//Left shock wave
		rho_starL=rho_L*(P_star/p_L+(gamma_L-1.)/(gamma_L+1.))/(P_star/p_L*(gamma_L-1.)/(gamma_L+1.)+1.);
	if(P_star<=p_R)//Right rarefaction wave
		rho_starR=rho_R*pow(P_star/p_R,1./gamma_R);
	else//Right shock wave
		rho_starR=rho_R*(P_star/p_R+(gamma_R-1.)/(gamma_R+1.))/(P_star/p_R*(gamma_R-1.)/(gamma_R+1.)+1.);

  return gap;
}
#endif
