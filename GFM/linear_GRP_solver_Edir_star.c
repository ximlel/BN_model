#include <math.h>
#include <stdio.h>
#include "mex.h"


double Riemann_solver_exact(double * U_star, double * P_star, double gammaL, double gammaR, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double eps, double tol, int N)
{
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
  if(p_R > p_L) // (u_LR,p_R) lies on the shock branch of I1
  {
    delta_p = p_R - p_L;
    u_LR = sqrt(1.0 + nuL*delta_p/p_L);
    u_LR = delta_p * c_L / gammaL / p_L / u_LR;
    u_LR = u_L - u_LR;
  }
  else // (u_LR,p_R) lies on the rarefaction branch of I1
  {
    u_LR = pow(p_R/p_L, muL) - 1.0;
    u_LR = 2.0 * c_L * u_LR / (gammaL-1.0);
    u_LR = u_L - u_LR;
  }
  //find out where (u_RL,p_L) lies on the curve of RIGHT state
  if(p_L > p_R) // (u_RL, p_L) lies on the shock branch of I3
  {
    delta_p = p_L - p_R;
    u_RL = sqrt(1.0 + nuR*delta_p/p_R);
    u_RL = delta_p * c_R / gammaR / p_R / u_RL;
    u_RL = u_R + u_RL;
  }
  else // (u_RL, p_L) lies on the rarefaction branch of I3
  {
    u_RL = pow(p_L/p_R, muR) - 1.0;
    u_RL = 2.0 * c_R * u_RL / (gammaR-1.0);
    u_RL = u_R + u_RL;
  }
  if(u_LR > u_R+eps)
    CRW[1] = 0;
  else
    CRW[1] = 1;
  if(u_RL > u_L-eps)
    CRW[0] = 1;
  else
    CRW[0] = 0;

  //======one step of the Newton ietration to get the intersection point of I1 and I3====
  k1 = -c_L / p_L / gammaL;//the (p,u)-tangent slope on I1 at (u_L,p_L), i.e. [du/dp](p_L)
  k3 =  c_R / p_R / gammaR;//the (p,u)-tangent slope on I3 at (u_R,p_R), i.e. [du/dp](p_R)
  //the intersect of (u-u_L)=k1*(p-p_L) and (u-u_R)=k3*(p-p_R)
  p_INT = (k1*p_L - k3*p_R - u_L + u_R) / (k1 - k3);
  if(p_INT < 0)
    p_INT = (p_L<p_R)? p_L : p_R;
  p_INT = 0.5*p_INT;

  //=======compute the gap between U^n_R and U^n_L(see Appendix C)=======
  if(p_INT > p_L)
  {
    delta_p = p_INT - p_L;
    v_L = sqrt(1.0 + nuL*delta_p/p_L);
    v_L = delta_p * c_L / gammaL / p_L / v_L;
    v_L = u_L - v_L;
  }
  else
  {
    v_L = pow(p_INT/p_L, muL) - 1.0;
    v_L = 2.0 * c_L * v_L / (gammaL-1.0);
    v_L = u_L - v_L;
  }
  if(p_INT > p_R)
  {
    delta_p = p_INT - p_R;
    v_R = sqrt(1.0 + nuR*delta_p/p_R);
    v_R = delta_p * c_R / gammaR / p_R / v_R;
    v_R = u_R + v_R;
  }
  else
  {
    //dbg = pow(p_INT/p_R, mu);
    v_R = pow(p_INT/p_R, muR) - 1.0;
    v_R = 2.0 * c_R * v_R / (gammaR-1.0);
    v_R = u_R + v_R;
  }
  gap = fabs(v_L - v_R);


  //=======THE NEWTON ITERATION=======
  while((gap > tol) && (n != N))
  {
    //the (p,u)-tangent slope on I1 at (v_L,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_L)
    {
      delta_p = p_INT - p_L;
      temp1 = 1.0 / sqrt(1.0 + nuL*delta_p/p_L);
      temp2 = c_L / gammaL / p_L;
      temp3 = 0.5 * temp2 * nuL / p_L;
      k1 = temp3*delta_p*pow(temp1,3.0) - temp2*temp1;
    }
    else
    {
      temp2 = c_L / gammaL / p_L;
      temp1 = 1.0 / pow(p_INT/p_L, nuL);
      k1 = -temp1 * temp2;
    }
    //the (p,u)-tangent slope on I3 at (v_R,p_INT), i.e. [du/dp](p_INT)
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      temp1 = 1.0 / sqrt(1.0 + nuR*delta_p/p_R);
      temp2 = c_R / gammaR / p_R;
      temp3 = 0.5 * temp2 * nuR / p_R;
      k3 = temp2*temp1 - temp3*delta_p*pow(temp1,3.0);
    }
    else
    {
      temp2 = c_R / gammaR / p_R;
      temp1 = 1.0 / pow(p_INT/p_R, nuR);
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
      v_L = sqrt(1.0 + nuL*delta_p/p_L);
      v_L = delta_p * c_L / gammaL / p_L / v_L;
      v_L = u_L - v_L;
    }
    else
    {
      v_L = pow(p_INT/p_L, muL) - 1.0;
      v_L = 2.0 * c_L * v_L / (gammaL-1.0);
      v_L = u_L - v_L;
    }
    if(p_INT > p_R)
    {
      delta_p = p_INT - p_R;
      v_R = sqrt(1.0 + nuR*delta_p/p_R);
      v_R = delta_p * c_R / gammaR / p_R / v_R;
      v_R = u_R + v_R;
    }
    else
    {
      v_R = pow(p_INT/p_R, muR) - 1.0;
      v_R = 2.0 * c_R * v_R / (gammaR-1.0);
      v_R = u_R + v_R;
    }

    gap = fabs(v_L - v_R);
  }
//printf("n=%d\tN=%d\n", n, N);

  u_INT = k1*(v_R-v_L)/(k1-k3)+v_L;

  *P_star = p_INT;
  *U_star = u_INT;

  return gap;
}



/*
atc=1.0/0.0(inf) acoustic approximation
atc=eps          Q1D GRP solver(nonlinear + acoustic case)
atc=-0.0         Q1D GRP solver(only nonlinear case)
atc,d_,t_=-0.0   exact Riemann solver
atc=eps,t_=-0.0  P1D GRP solver
*/

void linear_GRP_solver_Edir_Q1D
(double *wave_speed, double *D, double *U, double *U_star, double *D_star, double *input)
{
 double lambda_u=0.0, lambda_v=0.0;
 double rho_L, rho_R, d_rho_L, d_rho_R, t_rho_L, t_rho_R, u_L, u_R, d_u_L, d_u_R, t_u_L, t_u_R, v_L, v_R, d_v_L, d_v_R, t_v_L, t_v_R, p_L, p_R, d_p_L, d_p_R, t_p_L, t_p_R, phi_L, phi_R, d_phi_L, d_phi_R, t_phi_L, t_phi_R;
 double z_L=0.0, z_R=0.0, d_z_L=0.0, d_z_R=0.0, t_z_L=0.0, t_z_R=0.0;
 double gammaL, gammaR, eps, atc;
 rho_L=input[0];
 rho_R=input[1];
 d_rho_L=input[2];
 d_rho_R=input[3];
 t_rho_L=input[4];
 t_rho_R=input[5];
 u_L=input[6];
 u_R=input[7];
 d_u_L=input[8];
 d_u_R=input[9];
 t_u_L=input[10];
 t_u_R=input[11];
 v_L=input[12];
 v_R=input[13];
 d_v_L=input[14];
 d_v_R=input[15];
 t_v_L=input[16];
 t_v_R=input[17];
 p_L=input[18];
 p_R=input[19];
 d_p_L=input[20];
 d_p_R=input[21];
 t_p_L=input[22];
 t_p_R=input[23];
 phi_L=input[24];
 phi_R=input[25];
 d_phi_L=input[26];
 d_phi_R=input[27];
 t_phi_L=input[28];
 t_phi_R=input[29];
 gammaL=input[30];
 gammaR=input[31];
 eps=input[32];
 atc=input[33];

	int CRW[2]={0};
	double dist;
	double c_L, c_R, C, c_frac;

	double d_Phi, d_Psi, TdS, VAR;
	double D_rho, D_u, D_v, D_p, D_z, D_phi, T_rho, T_u, T_v, T_p, T_z, T_phi; 
	double u_star, p_star, rho_star_L, rho_star_R, c_star_L, c_star_R;

	double H1, H2, H3;
	double a_L, b_L, d_L, a_R, b_R, d_R, detA;
	double L_u, L_p, L_rho;

	double u_t_mat, p_t_mat;
	double SmUs, SmUL, SmUR;
  
	const double zetaL = (gammaL-1.0)/(gammaL+1.0);
	const double zetaR = (gammaR-1.0)/(gammaR+1.0);
 
	double rho_x;
	double g_rho, g_u, g_p, f;

	double speed_L, speed_R;

	c_L = sqrt(gammaL * p_L / rho_L);
	c_R = sqrt(gammaR * p_R / rho_R);


	//=========non-acoustic case==========
	Riemann_solver_exact(&u_star, &p_star, gammaL, gammaR, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, eps, 500);

	if(CRW[0])
		{
			rho_star_L = rho_L*pow(p_star/p_L, 1.0/gammaL);
			c_star_L = c_L*pow(p_star/p_L, 0.5*(gammaL-1.0)/gammaL);
			speed_L = u_L - c_L;
		}
	else
		{
			rho_star_L = rho_L*(p_star+zetaL*p_L)/(p_L+zetaL*p_star);
			c_star_L = sqrt(gammaL * p_star / rho_star_L);
			speed_L = u_L - c_L*sqrt(0.5*((gammaL+1.0)*(p_star/p_L) + (gammaL-1.0))/gammaL);
		}
	if(CRW[1])
		{
			rho_star_R = rho_R*pow(p_star/p_R,1.0/gammaR);
			c_star_R = c_R*pow(p_star/p_R, 0.5*(gammaR-1.0)/gammaR);
			speed_R = u_R + c_R;
		}
	else
		{
			rho_star_R = rho_R*(p_star+zetaR*p_R)/(p_R+zetaR*p_star);
			c_star_R = sqrt(gammaR * p_star / rho_star_R);
			speed_R = u_R + c_R*sqrt(0.5*((gammaR+1.0)*(p_star/p_R) + (gammaR-1.0))/gammaR);
		}
	wave_speed[0] = speed_L;
	wave_speed[1] = speed_R;


				//--non-sonic case--
				{

					//determine a_L, b_L and d_L
					if(CRW[0]) //the 1-wave is a CRW
						{
							a_L = 1.0;
							b_L = 1.0 / rho_star_L / c_star_L;
							c_frac = c_star_L/c_L;
							TdS = (d_p_L - d_rho_L*c_L*c_L)/(gammaL-1.0)/rho_L;
							d_Psi = d_u_L + (gammaL*d_p_L/c_L - c_L*d_rho_L)/(gammaL-1.0)/rho_L;
							d_L = ((1.0+zetaL)*pow(c_frac, 0.5/zetaL) + zetaL*pow(c_frac, (1.0+zetaL)/zetaL));
							d_L = d_L/(1.0+2.0*zetaL) * TdS;
							d_L = d_L - c_L*pow(c_frac, 0.5/zetaL) * d_Psi;
						}
					else //the 1-wave is a shock
						{
							SmUs = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
							SmUL = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

							VAR = sqrt((1-zetaL)/(rho_L*(p_star+zetaL*p_L)));

							H1 =  0.5*VAR * (p_star+(1.0+2.0*zetaL)*p_L)/(p_star+zetaL*p_L);
							H2 = -0.5*VAR * ((2.0+zetaL)*p_star + zetaL*p_L)/(p_star+zetaL*p_L);
							H3 = -0.5*VAR * (p_star-p_L) / rho_L;

							L_p = -1.0/rho_L - SmUL*H2;
							L_u = SmUL + rho_L*(c_L*c_L*H2 + H3);
							L_rho = -SmUL * H3;

							a_L = 1.0 - rho_star_L* SmUs * H1;
							b_L = -SmUs/(rho_star_L*c_star_L*c_star_L)+ H1;
							d_L = L_rho*d_rho_L + L_u*d_u_L + L_p*d_p_L;
						}
					//determine a_R, b_R and d_R
					if(CRW[1]) //the 3-wave is a CRW
						{
							a_R = 1.0;
							b_R = -1.0 / rho_star_R / c_star_R;
							c_frac = c_star_R/c_R;
							TdS = (d_p_R - d_rho_R*c_R*c_R)/(gammaR-1.0)/rho_R;
							d_Phi = d_u_R - (gammaR*d_p_R/c_R - c_R*d_rho_R)/(gammaR-1.0)/rho_R;
							d_R = ((1.0+zetaR)*pow(c_frac, 0.5/zetaR) + zetaR*pow(c_frac, (1.0+zetaR)/zetaR));
							d_R = d_R/(1.0+2.0*zetaR) * TdS;
							d_R = d_R + c_R*pow(c_frac, 0.5/zetaR)*d_Phi;
						}
					else //the 3-wave is a shock
						{
							SmUs = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
							SmUR = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

							VAR  = sqrt((1.0-zetaR)/(rho_R*(p_star+zetaR*p_R)));

							H1 = 0.5* VAR * (p_star+(1+2.0*zetaR)*p_R)/(p_star+zetaR*p_R);
							H2 = -0.5*VAR * ((2.0+zetaR)*p_star+zetaR*p_R)/(p_star+zetaR*p_R);
							H3 = -0.5*(p_star-p_R)* VAR /rho_R;

							L_p = -1.0/rho_R + SmUR*H2;
							L_u = SmUR - rho_R*(c_R*c_R*H2 + H3);
							L_rho = SmUR * H3;

							a_R = 1.0 +rho_star_R* SmUs * H1;
							b_R = -(SmUs/(rho_star_R*c_star_R*c_star_R) + H1);
							d_R = L_rho*d_rho_R + L_u*d_u_R + L_p*d_p_R;
						}

					detA = a_L*b_R - b_L*a_R;
					u_t_mat = (b_R*d_L - b_L*d_R)/detA;
					p_t_mat = (a_L*d_R - a_R*d_L)/detA;


//					if(u_star < lambda_u) //the direction is between the contact discontinuety and the 3-wave
						{
							U[0] = rho_star_R;
							U[1] =   u_star;
							U[2] =   v_R;
							U[3] =   p_star;
							U[4] =   z_R;
							U[5] = phi_R;
							C = c_star_R;
					//already total D!
					D[1] = u_t_mat + (u_star-lambda_u)/U[0]/C/C * p_t_mat;
					D[3] = p_t_mat + (u_star-lambda_u)*U[0] * u_t_mat;
					D_star[4]=-p_t_mat/U[0]/C/C;
					D_star[5]=-U[0] * u_t_mat;
							if(CRW[1]) //the 3-wave is a CRW
								{
									//already total D!
									D[0] = rho_star_R*(u_star-lambda_u)*pow(c_star_R/c_R, (1.0+zetaR)/zetaR)*(d_p_R - d_rho_R*c_R*c_R)/rho_R;
									D[0] = (D[0] + D[3]) / c_star_R/c_star_R;

									D[2] = -U[1]*d_v_R*U[0]/rho_R;
									D[2] = D[2] + lambda_u*d_v_R;
									D[4] = -U[1]*d_z_R*U[0]/rho_R;
									D[4] = D[4] + lambda_u*d_z_R;
									D[5] = -U[1]*d_phi_R*U[0]/rho_R;
									D[5] = D[5] + lambda_u*d_phi_R;
								}
							else //the 3-wave is a shock
								{
									SmUs = sqrt(0.5*((gammaR+1.0)*p_R   + (gammaR-1.0)*p_star)/rho_star_R);
									SmUR = sqrt(0.5*((gammaR+1.0)*p_star+ (gammaR-1.0)*p_R   )/rho_R);

									VAR = p_R + zetaR*p_star;
									H1 = rho_R * p_R    * (1.0 - zetaR*zetaR) / VAR/VAR;
									H2 = rho_R * p_star * (zetaR*zetaR - 1.0) / VAR/VAR;
									H3 = (p_star + zetaR*p_R)/VAR;

									L_rho = SmUR * H3 * d_rho_R;
									L_u = -rho_R * (H2*c_R*c_R + H3) * d_u_R;
									L_p = H2 * SmUR * d_p_R;

									D[0] = ((u_star+SmUs)/c_star_R/c_star_R - u_star*H1)*p_t_mat + rho_star_R*u_star*SmUs*H1*u_t_mat;
									D[0] = (D[0] - u_star*(L_p+L_rho+L_u)) / SmUs;

									f = SmUR*(H2*d_p_R + H3*d_rho_R) - rho_R*(H2*c_R*c_R+H3)*d_u_R;
									rho_x = (f + H1*(p_t_mat - rho_star_R*SmUs*u_t_mat) - D[0]) / (SmUR+u_R);//shk_spd;
									D[0] = D[0] + lambda_u*rho_x;
									D_star[2]=rho_x;

									D[2] = -U[1] * SmUR * d_v_R / SmUs;
									D[2] = D[2] + lambda_u*d_v_R;
									D[4] = -U[1] * SmUR * d_z_R / SmUs;
									D[4] = D[4] + lambda_u*d_z_R;
									D[5] = -U[1] * SmUR * d_phi_R / SmUs;
									D[5] = D[5] + lambda_u*d_phi_R;
								}
						}
//					else //the direction is between the 1-wave and the contact discontinuety
						{
							U[0] = rho_star_L;
							U[1] =   u_star;
							U[2] =   v_L;
							U[3] =   p_star;
							U[4] =   z_L;
							U[5] = phi_L;
							C = c_star_L;
					//already total D!
					D[1] = u_t_mat + (u_star-lambda_u)/U[0]/C/C * p_t_mat;
					D[3] = p_t_mat + (u_star-lambda_u)*U[0] * u_t_mat;
					D_star[1]=-p_t_mat/U[0]/C/C;
					D_star[3]=-U[0] * u_t_mat;
							if(CRW[0]) //the 1-wave is a CRW
								{
									//already total D!
									D[0] = rho_star_L*(u_star-lambda_u)*pow(c_star_L/c_L, (1.0+zetaL)/zetaL)*(d_p_L - d_rho_L*c_L*c_L)/rho_L;
									D[0] = (D[0] + D[3]) / c_star_L/c_star_L;

									D[2] = -U[1]*d_v_L*U[0]/rho_L;
									D[2] = D[2] + lambda_u*d_v_L;
									D[4] = -U[1]*d_z_L*U[0]/rho_L;
									D[4] = D[4] + lambda_u*d_z_L;
									D[5] = -U[1]*d_phi_L*U[0]/rho_L;
									D[5] = D[5] + lambda_u*d_phi_L;
								}
							else //the 1-wave is a shock
								{
									SmUs = -sqrt(0.5*((gammaL+1.0)*p_L   +(gammaL-1.0)*p_star)/rho_star_L);
									SmUL = -sqrt(0.5*((gammaL+1.0)*p_star+(gammaL-1.0)*p_L   )/rho_L);

									VAR = p_L + zetaL*p_star;

									H1 = rho_L * p_L    * (1.0 - zetaL*zetaL) / VAR/VAR;
									H2 = rho_L * p_star * (zetaL*zetaL - 1.0) / VAR/VAR;
									H3 = (p_star + zetaL*p_L)/VAR;

									L_rho = SmUL * H3 * d_rho_L;
									L_u = -rho_L*(H2*c_L*c_L + H3) * d_u_L;
									L_p = H2 * SmUL * d_p_L;

									D[0] = ((u_star+SmUs)/c_star_L/c_star_L - H1*u_star)*p_t_mat + rho_star_L*u_star*SmUs*H1*u_t_mat;
									D[0] = (D[0] - u_star*(L_p+L_rho+L_u))/ SmUs;

									f = SmUL*(H2*d_p_L + H3*d_rho_L) - rho_L*(H2*c_L*c_L+H3)*d_u_L;
									rho_x = (f + H1*(p_t_mat - rho_star_L*SmUs*u_t_mat) - D[0]) / (SmUL+u_L);
									D[0] = D[0] + lambda_u*rho_x;
									D_star[0]=rho_x;

									D[2] = -U[1] * SmUL * d_v_L / SmUs;
									D[2] = D[2] + lambda_u*d_v_L;
									D[4] = -U[1] * SmUL * d_z_L / SmUs;
									D[4] = D[4] + lambda_u*d_z_L;
									D[5] = -U[1] * SmUL * d_phi_L / SmUs;
									D[5] = D[5] + lambda_u*d_phi_L;
								}
						}
					//--end of non-sonic case--
				}
			T_p = 0.5*((t_u_L*(U[0]*C) + t_p_L) - (t_u_R*(U[0]*C) - t_p_R));
			T_u = 0.5*(t_u_L + t_p_L/(U[0]*C) + t_u_R - t_p_R/(U[0]*C));
			if (u_star > lambda_u)
				{
					T_rho = t_rho_L - t_p_L/(C*C) + T_p/(C*C);
					D[0] = D[0] - (U[2]-lambda_v)*T_rho - U[0]*t_v_L;
					D[1] = D[1] - (U[2]-lambda_v)*T_u;
					D[2] = D[2] - (U[2]-lambda_v)*t_v_L - T_p/U[0];
					D[3] = D[3] - (U[2]-lambda_v)*T_p   - U[0]*C*C*t_v_L;
					D[4] = D[4] - (U[2]-lambda_v)*t_z_L;
					D[5] = D[5] - (U[2]-lambda_v)*t_phi_L;							
				}
			else
				{
					T_rho = t_rho_R - t_p_R/(C*C) + T_p/(C*C);
					D[0] = D[0] - (U[2]-lambda_v)*T_rho - U[0]*t_v_R;
					D[1] = D[1] - (U[2]-lambda_v)*T_u;
					D[2] = D[2] - (U[2]-lambda_v)*t_v_R - T_p/U[0];
					D[3] = D[3] - (U[2]-lambda_v)*T_p   - U[0]*C*C*t_v_R;
					D[4] = D[4] - (U[2]-lambda_v)*t_z_R;
					D[5] = D[5] - (U[2]-lambda_v)*t_phi_R;
				}

	U_star[0] = rho_star_L;
	U_star[1] = u_star;
	U_star[2] = rho_star_R;
	U_star[3] = p_star;
	U_star[4] = c_star_L;
	U_star[5] = c_star_R;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    double *M;
    int m,n;
    M = mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
     
    plhs[0] = mxCreateDoubleMatrix(26,1,mxREAL);
    double *A;
    A = mxGetPr(plhs[0]);
     
    linear_GRP_solver_Edir_Q1D(A,A+2,A+8,A+14,A+20,M);
}
