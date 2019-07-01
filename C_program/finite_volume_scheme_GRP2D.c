#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
#include "../include/finite_volume.h"
#include "../include/Riemann_solver.h"
#include "../include/file_io.h"

#define pi (4.*atan(1.0))

static double minmod3(const double a, double b, const double c)
{
	const double Alpha=1.5;
	if (a>0.0&&b>0.0&&c>0.0)
		return fmin(fmin(Alpha*a,Alpha*c),b);
	if (a<0.0&&b<0.0&&c<0.0)
		return fmax(fmax(Alpha*a,Alpha*c),b);
	else
		return 0.0;
}

static double minmod2(const double a, double b)
{
	if (a>0.0&&b>0.0)
		return fmin(a,b);
	if (a<0.0&&b<0.0)
		return fmax(a,b);
	else
		return 0.0;
}

static void slope_limiter_GRP(int n_x, double Z_gC[][n_x], double Z_I_gxL[][n_x+1], double Z_I_gxR[][n_x+1], double Z_I_gyL[][n_x+1], double Z_I_gyR[][n_x+1], double Z_gx[][n_x], double Z_gy[][n_x],
							  double RHO_gC[][n_x], double P_gC[][n_x], double U_gC[][n_x], double V_gC[][n_x], double RHO_I_gxL[][n_x+1], double RHO_I_gxR[][n_x+1], double P_I_gxL[][n_x+1], double P_I_gxR[][n_x+1], double U_I_gxL[][n_x+1], double U_I_gxR[][n_x+1], double V_I_gxL[][n_x+1], double V_I_gxR[][n_x+1], double RHO_I_gyL[][n_x+1], double RHO_I_gyR[][n_x+1], double P_I_gyL[][n_x+1], double P_I_gyR[][n_x+1], double U_I_gyL[][n_x+1], double U_I_gyR[][n_x+1], double V_I_gyL[][n_x+1], double V_I_gyR[][n_x+1], double RHO_gx[][n_x], double P_gx[][n_x], double U_gx[][n_x], double V_gx[][n_x], double RHO_gy[][n_x], double P_gy[][n_x], double U_gy[][n_x], double V_gy[][n_x],
							  double RHO_lC[][n_x], double P_lC[][n_x], double U_lC[][n_x], double V_lC[][n_x], double RHO_I_lxL[][n_x+1], double RHO_I_lxR[][n_x+1], double P_I_lxL[][n_x+1], double P_I_lxR[][n_x+1], double U_I_lxL[][n_x+1], double U_I_lxR[][n_x+1], double V_I_lxL[][n_x+1], double V_I_lxR[][n_x+1], double RHO_I_lyL[][n_x+1], double RHO_I_lyR[][n_x+1], double P_I_lyL[][n_x+1], double P_I_lyR[][n_x+1], double U_I_lyL[][n_x+1], double U_I_lyR[][n_x+1], double V_I_lyL[][n_x+1], double V_I_lyR[][n_x+1], double RHO_lx[][n_x], double P_lx[][n_x], double U_lx[][n_x], double V_lx[][n_x], double RHO_ly[][n_x], double P_ly[][n_x], double U_ly[][n_x], double V_ly[][n_x])
{
	const int n_y = (int)config[14];
	for(int i = 0; i < n_y; ++i)
		for(int j = 0; j < n_x; ++j)
			{
				Z_gx[i][j]   = (Z_I_gxL[i][j+1]   - Z_I_gxR[i][j])  /config[10];
				Z_gy[i][j]   = (Z_I_gyL[i+1][j]   - Z_I_gyR[i][j])  /config[11];
				RHO_gx[i][j] = (RHO_I_gxL[i][j+1] - RHO_I_gxR[i][j])/config[10];
				RHO_gy[i][j] = (RHO_I_gyL[i+1][j] - RHO_I_gyR[i][j])/config[11];
				U_gx[i][j]   = (U_I_gxL[i][j+1]   - U_I_gxR[i][j])  /config[10];
				U_gy[i][j]   = (U_I_gyL[i+1][j]   - U_I_gyR[i][j])  /config[11];
				V_gx[i][j]   = (V_I_gxL[i][j+1]   - V_I_gxR[i][j])  /config[10];
				V_gy[i][j]   = (V_I_gyL[i+1][j]   - V_I_gyR[i][j])  /config[11];
				P_gx[i][j]   = (P_I_gxL[i][j+1]   - P_I_gxR[i][j])  /config[10];
				P_gy[i][j]   = (P_I_gyL[i+1][j]   - P_I_gyR[i][j])  /config[11];
				RHO_lx[i][j] = (RHO_I_lxL[i][j+1] - RHO_I_lxR[i][j])/config[10];
				RHO_ly[i][j] = (RHO_I_lyL[i+1][j] - RHO_I_lyR[i][j])/config[11];
				U_lx[i][j]   = (U_I_lxL[i][j+1]   - U_I_lxR[i][j])  /config[10];
				U_ly[i][j]   = (U_I_lyL[i+1][j]   - U_I_lyR[i][j])  /config[11];
				V_lx[i][j]   = (V_I_lxL[i][j+1]   - V_I_lxR[i][j])  /config[10];
				V_ly[i][j]   = (V_I_lyL[i+1][j]   - V_I_lyR[i][j])  /config[11];
				P_lx[i][j]   = (P_I_lxL[i][j+1]   - P_I_lxR[i][j])  /config[10];
				P_ly[i][j]   = (P_I_lyL[i+1][j]   - P_I_lyR[i][j])  /config[11];

				int iL,iR,jL,jR;
				if (i==0)
					{
						iL=n_y-1;
						iR=i+1;
					}
				else if (i==n_y-1)
					{
						iL=i-1;
						iR=0;
					}
				else
					{
						iL=i-1;
						iR=i+1;
					}
				if (j==0)
					{
						jL=j;
						jR=j+1;
					}
				else if (j==n_x-1)
					{
						jL=j-1;
						jR=j;
					}
				else
					{
						jL=j-1;
						jR=j+1;
					}

				Z_gx[i][j] = minmod3((Z_gC[i][j]-Z_gC[i][jL])/config[10],Z_gx[i][j],(Z_gC[i][jR]-Z_gC[i][j])/config[10]);
				Z_gy[i][j] = minmod3((Z_gC[i][j]-Z_gC[iL][j])/config[11],Z_gy[i][j],(Z_gC[iR][j]-Z_gC[i][j])/config[11]);
				RHO_gx[i][j] = minmod3((RHO_gC[i][j]-RHO_gC[i][jL])/config[10],RHO_gx[i][j],(RHO_gC[i][jR]-RHO_gC[i][j])/config[10]);
				RHO_gy[i][j] = minmod3((RHO_gC[i][j]-RHO_gC[iL][j])/config[11],RHO_gy[i][j],(RHO_gC[iR][j]-RHO_gC[i][j])/config[11]);
				U_gx[i][j] = minmod3((U_gC[i][j]-U_gC[i][jL])/config[10],U_gx[i][j],(U_gC[i][jR]-U_gC[i][j])/config[10]);
				U_gy[i][j] = minmod3((U_gC[i][j]-U_gC[iL][j])/config[11],U_gy[i][j],(U_gC[iR][j]-U_gC[i][j])/config[11]);
				V_gx[i][j] = minmod3((V_gC[i][j]-V_gC[i][jL])/config[10],V_gx[i][j],(V_gC[i][jR]-V_gC[i][j])/config[10]);
				V_gy[i][j] = minmod3((V_gC[i][j]-V_gC[iL][j])/config[11],V_gy[i][j],(V_gC[iR][j]-V_gC[i][j])/config[11]);
				P_gx[i][j] = minmod3((P_gC[i][j]-P_gC[i][jL])/config[10],P_gx[i][j],(P_gC[i][jR]-P_gC[i][j])/config[10]);
				P_gy[i][j] = minmod3((P_gC[i][j]-P_gC[iL][j])/config[11],P_gy[i][j],(P_gC[iR][j]-P_gC[i][j])/config[11]);
				RHO_lx[i][j] = minmod3((RHO_lC[i][j]-RHO_lC[i][jL])/config[10],RHO_lx[i][j],(RHO_lC[i][jR]-RHO_lC[i][j])/config[10]);
				RHO_ly[i][j] = minmod3((RHO_lC[i][j]-RHO_lC[iL][j])/config[11],RHO_ly[i][j],(RHO_lC[iR][j]-RHO_lC[i][j])/config[11]);
				U_lx[i][j] = minmod3((U_lC[i][j]-U_lC[i][jL])/config[10],U_lx[i][j],(U_lC[i][jR]-U_lC[i][j])/config[10]);
				U_ly[i][j] = minmod3((U_lC[i][j]-U_lC[iL][j])/config[11],U_ly[i][j],(U_lC[iR][j]-U_lC[i][j])/config[11]);
				V_lx[i][j] = minmod3((V_lC[i][j]-V_lC[i][jL])/config[10],V_lx[i][j],(V_lC[i][jR]-V_lC[i][j])/config[10]);
				V_ly[i][j] = minmod3((V_lC[i][j]-V_lC[iL][j])/config[11],V_ly[i][j],(V_lC[iR][j]-V_lC[i][j])/config[11]);
				P_lx[i][j] = minmod3((P_lC[i][j]-P_lC[i][jL])/config[10],P_lx[i][j],(P_lC[i][jR]-P_lC[i][j])/config[10]);
				P_ly[i][j] = minmod3((P_lC[i][j]-P_lC[iL][j])/config[11],P_ly[i][j],(P_lC[iR][j]-P_lC[i][j])/config[11]);
			}
}

static void slope_limiter3_GRP(int n_x, double Z_gC[][n_x], double Z_I_gxL[][n_x+1], double Z_I_gxR[][n_x+1], double Z_I_gyL[][n_x+1], double Z_I_gyR[][n_x+1], double Z_gx[][n_x], double Z_gy[][n_x],
							   double RHO_gC[][n_x], double P_gC[][n_x], double U_gC[][n_x], double V_gC[][n_x], double RHO_I_gxL[][n_x+1], double RHO_I_gxR[][n_x+1], double P_I_gxL[][n_x+1], double P_I_gxR[][n_x+1], double U_I_gxL[][n_x+1], double U_I_gxR[][n_x+1], double V_I_gxL[][n_x+1], double V_I_gxR[][n_x+1], double RHO_I_gyL[][n_x+1], double RHO_I_gyR[][n_x+1], double P_I_gyL[][n_x+1], double P_I_gyR[][n_x+1], double U_I_gyL[][n_x+1], double U_I_gyR[][n_x+1], double V_I_gyL[][n_x+1], double V_I_gyR[][n_x+1], double RHO_gx[][n_x], double P_gx[][n_x], double U_gx[][n_x], double V_gx[][n_x], double RHO_gy[][n_x], double P_gy[][n_x], double U_gy[][n_x], double V_gy[][n_x],
							   double RHO_lC[][n_x], double P_lC[][n_x], double U_lC[][n_x], double V_lC[][n_x], double RHO_I_lxL[][n_x+1], double RHO_I_lxR[][n_x+1], double P_I_lxL[][n_x+1], double P_I_lxR[][n_x+1], double U_I_lxL[][n_x+1], double U_I_lxR[][n_x+1], double V_I_lxL[][n_x+1], double V_I_lxR[][n_x+1], double RHO_I_lyL[][n_x+1], double RHO_I_lyR[][n_x+1], double P_I_lyL[][n_x+1], double P_I_lyR[][n_x+1], double U_I_lyL[][n_x+1], double U_I_lyR[][n_x+1], double V_I_lyL[][n_x+1], double V_I_lyR[][n_x+1], double RHO_lx[][n_x], double P_lx[][n_x], double U_lx[][n_x], double V_lx[][n_x], double RHO_ly[][n_x], double P_ly[][n_x], double U_ly[][n_x], double V_ly[][n_x])
{
	const int n_y = (int)config[14];
	int iL,iR,jL,jR;
	for(int i = 0; i < n_y; ++i)
		for(int j = 0; j < n_x; ++j)
			{
				if (i==0)
					{
						iL=n_y-1;
						iR=i+1;
					}
				else if (i==n_y-1)
					{
						iL=i-1;
						iR=0;
					}
				else
					{
						iL=i-1;
						iR=i+1;
					}
				if (j==0)
					{
						jL=j;
						jR=j+1;
					}
				else if (j==n_x-1)
					{
						jL=j-1;
						jR=j;
					}
				else
					{
						jL=j-1;
						jR=j+1;
					}

				Z_gx[i][j] = (Z_gC[i][jR]-Z_gC[i][jL])/config[10]/2;
				Z_gy[i][j] = (Z_gC[iR][j]-Z_gC[iL][j])/config[11]/2;
				RHO_gx[i][j] = (RHO_gC[i][jR]-RHO_gC[i][jL])/config[10]/2;
				RHO_gy[i][j] = (RHO_gC[iR][j]-RHO_gC[iL][j])/config[11]/2;
				U_gx[i][j] = (U_gC[i][jR]-U_gC[i][jL])/config[10]/2;
				U_gy[i][j] = (U_gC[iR][j]-U_gC[iL][j])/config[11]/2;
				V_gx[i][j] = (V_gC[i][jR]-V_gC[i][jL])/config[10]/2;
				V_gy[i][j] = (V_gC[iR][j]-V_gC[iL][j])/config[11]/2;
				P_gx[i][j] = (P_gC[i][jR]-P_gC[i][jL])/config[10]/2;
				P_gy[i][j] = (P_gC[iR][j]-P_gC[iL][j])/config[11]/2;
				RHO_lx[i][j] = (RHO_lC[i][jR]-RHO_lC[i][jL])/config[10]/2;
				RHO_ly[i][j] = (RHO_lC[iR][j]-RHO_lC[iL][j])/config[11]/2;
				U_lx[i][j] = (U_lC[i][jR]-U_lC[i][jL])/config[10]/2;
				U_ly[i][j] = (U_lC[iR][j]-U_lC[iL][j])/config[11]/2;
				V_lx[i][j] = (V_lC[i][jR]-V_lC[i][jL])/config[10]/2;
				V_ly[i][j] = (V_lC[iR][j]-V_lC[iL][j])/config[11]/2;
				P_lx[i][j] = (P_lC[i][jR]-P_lC[i][jL])/config[10]/2;
				P_ly[i][j] = (P_lC[iR][j]-P_lC[iL][j])/config[11]/2;

				Z_gx[i][j] = minmod2((Z_gC[i][jR]-Z_gC[i][j])/config[10],(Z_gC[i][j]-Z_gC[i][jL])/config[10]);
				Z_gy[i][j] = minmod2((Z_gC[iR][j]-Z_gC[i][j])/config[11],(Z_gC[i][j]-Z_gC[iL][j])/config[11]);
				RHO_gx[i][j] = minmod2((RHO_gC[i][jR]-RHO_gC[i][j])/config[10],(RHO_gC[i][j]-RHO_gC[i][jL])/config[10]);
				RHO_gy[i][j] = minmod2((RHO_gC[iR][j]-RHO_gC[i][j])/config[11],(RHO_gC[i][j]-RHO_gC[iL][j])/config[11]);
				U_gx[i][j] = minmod2((U_gC[i][jR]-U_gC[i][j])/config[10],(U_gC[i][j]-U_gC[i][jL])/config[10]);
				U_gy[i][j] = minmod2((U_gC[iR][j]-U_gC[i][j])/config[11],(U_gC[i][j]-U_gC[iL][j])/config[11]);
				V_gx[i][j] = minmod2((V_gC[i][jR]-V_gC[i][j])/config[10],(V_gC[i][j]-V_gC[i][jL])/config[10]);
				V_gy[i][j] = minmod2((V_gC[iR][j]-V_gC[i][j])/config[11],(V_gC[i][j]-V_gC[iL][j])/config[11]);
				P_gx[i][j] = minmod2((P_gC[i][jR]-P_gC[i][j])/config[10],(P_gC[i][j]-P_gC[i][jL])/config[10]);
				P_gy[i][j] = minmod2((P_gC[iR][j]-P_gC[i][j])/config[11],(P_gC[i][j]-P_gC[iL][j])/config[11]);
				RHO_lx[i][j] = minmod2((RHO_lC[i][jR]-RHO_lC[i][j])/config[10],(RHO_lC[i][j]-RHO_lC[i][jL])/config[10]);
				RHO_ly[i][j] = minmod2((RHO_lC[iR][j]-RHO_lC[i][j])/config[11],(RHO_lC[i][j]-RHO_lC[iL][j])/config[11]);
				U_lx[i][j] = minmod2((U_lC[i][jR]-U_lC[i][j])/config[10],(U_lC[i][j]-U_lC[i][jL])/config[10]);
				U_ly[i][j] = minmod2((U_lC[iR][j]-U_lC[i][j])/config[11],(U_lC[i][j]-U_lC[iL][j])/config[11]);
				V_lx[i][j] = minmod2((V_lC[i][jR]-V_lC[i][j])/config[10],(V_lC[i][j]-V_lC[i][jL])/config[10]);
				V_ly[i][j] = minmod2((V_lC[iR][j]-V_lC[i][j])/config[11],(V_lC[i][j]-V_lC[iL][j])/config[11]);
				P_lx[i][j] = minmod2((P_lC[i][jR]-P_lC[i][j])/config[10],(P_lC[i][j]-P_lC[i][jL])/config[10]);
				P_ly[i][j] = minmod2((P_lC[iR][j]-P_lC[i][j])/config[11],(P_lC[i][j]-P_lC[iL][j])/config[11]);
			}
}

void finite_volume_scheme_GRP2D(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem)
{
	clock_t start_clock;
	double cpu_time = 0.0;

	//config[53] = 1.0;//R-K

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const double eps = config[4];
	struct cell_var cv = cell_mem_init(mv, FV);

	vol_comp(&cv, mv);

	cell_rel(&cv, mv);

	if (order > 1)
		cell_centroid(&cv, mv);

	printf("Grid has been constructed.\n");

	double tau; // the length of the time step
	double t_all = 0.0;

	const double delta_plot_t = 0.01;
	double plot_t = 0.0;
	int j, stop_step = 0, stop_t = 0;
	double eps_big = 1e-10, eps_big2 = 1e-3, eps3=1e-10;

	const double gamma_g=config[6], gamma_l=config[106];
	const int n_y = (int)config[14], n_x = (int)config[13];
	double Z_gC[n_y][n_x];
	double ZRHO_gC[n_y][n_x], RHO_gC[n_y][n_x], P_gC[n_y][n_x], U_gC[n_y][n_x], V_gC[n_y][n_x];
	double RHO_U_gC[n_y][n_x], RHO_V_gC[n_y][n_x], E_gC[n_y][n_x];
	double ZRHO_lC[n_y][n_x], RHO_lC[n_y][n_x], P_lC[n_y][n_x], U_lC[n_y][n_x], V_lC[n_y][n_x];
	double RHO_U_lC[n_y][n_x], RHO_V_lC[n_y][n_x], E_lC[n_y][n_x];
	double Z_g_2D[4];
	double U_RHO_g_2D[4], U_U_g_2D[4], U_V_g_2D[4], U_E_g_2D[4];
	double U_RHO_l_2D[4], U_U_l_2D[4], U_V_l_2D[4], U_E_l_2D[4];
	init_data_1(Z_g_2D, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D);
	for(int i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				if (j >=n_x/2 && i < n_y/2)
					Z_gC[i][j]=Z_g_2D[0];
				else if (j < n_x/2 && i < n_y/2)
					Z_gC[i][j]=Z_g_2D[1];
				else if (j < n_x/2 && i >=n_y/2)
					Z_gC[i][j]=Z_g_2D[2];
				else if (j >=n_x/2 && i >=n_y/2)
					Z_gC[i][j]=Z_g_2D[3];
			}
	for(int i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				if (j ==n_x/2 && i ==n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 0,1,2,3);
				else if (j < n_x/2 && i ==n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 2,2,3,3);
				else if (j > n_x/2 && i ==n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 1,1,4,4);
				else if (j ==n_x/2 && i < n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 1,1,2,2);
				else if (j ==n_x/2 && i > n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 3,3,4,4);				
				else if (j > n_x/2 && i < n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 0,0,0,0);
				else if (j < n_x/2 && i < n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 1,1,1,1);
				else if (j < n_x/2 && i > n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 2,2,2,2);
				else if (j > n_x/2 && i > n_y/2)
					U_init(n_x, ZRHO_gC,RHO_U_gC,RHO_V_gC,E_gC,ZRHO_lC,RHO_U_lC,RHO_V_lC,E_lC,i,j, U_RHO_g_2D, U_U_g_2D, U_V_g_2D, U_E_g_2D, U_RHO_l_2D, U_U_l_2D, U_V_l_2D, U_E_l_2D, 3,3,3,3);
			}
	double Z_gx[n_y][n_x], Z_gy[n_y][n_x], Z_I_gxL[n_y+1][n_x+1], Z_I_gxR[n_y+1][n_x+1], Z_I_gyL[n_y+1][n_x+1], Z_I_gyR[n_y+1][n_x+1];
	double RHO_gx[n_y][n_x], P_gx[n_y][n_x], U_gx[n_y][n_x], V_gx[n_y][n_x];
	double RHO_gy[n_y][n_x], P_gy[n_y][n_x], U_gy[n_y][n_x], V_gy[n_y][n_x];
	double RHO_I_gxL[n_y+1][n_x+1], P_I_gxL[n_y+1][n_x+1], U_I_gxL[n_y+1][n_x+1], V_I_gxL[n_y+1][n_x+1];
	double RHO_I_gxR[n_y+1][n_x+1], P_I_gxR[n_y+1][n_x+1], U_I_gxR[n_y+1][n_x+1], V_I_gxR[n_y+1][n_x+1];
	double RHO_I_gyL[n_y+1][n_x+1], P_I_gyL[n_y+1][n_x+1], U_I_gyL[n_y+1][n_x+1], V_I_gyL[n_y+1][n_x+1];
	double RHO_I_gyR[n_y+1][n_x+1], P_I_gyR[n_y+1][n_x+1], U_I_gyR[n_y+1][n_x+1], V_I_gyR[n_y+1][n_x+1];
	double RHO_F_gx[n_y+1][n_x+1], E_F_gx[n_y+1][n_x+1], U_F_gx[n_y+1][n_x+1], V_F_gx[n_y+1][n_x+1];
	double RHO_F_gy[n_y+1][n_x+1], E_F_gy[n_y+1][n_x+1], U_F_gy[n_y+1][n_x+1], V_F_gy[n_y+1][n_x+1];
	double RHO_lx[n_y][n_x], P_lx[n_y][n_x], U_lx[n_y][n_x], V_lx[n_y][n_x];
	double RHO_ly[n_y][n_x], P_ly[n_y][n_x], U_ly[n_y][n_x], V_ly[n_y][n_x];
	double RHO_I_lxL[n_y+1][n_x+1], P_I_lxL[n_y+1][n_x+1], U_I_lxL[n_y+1][n_x+1], V_I_lxL[n_y+1][n_x+1];
	double RHO_I_lxR[n_y+1][n_x+1], P_I_lxR[n_y+1][n_x+1], U_I_lxR[n_y+1][n_x+1], V_I_lxR[n_y+1][n_x+1];
	double RHO_I_lyL[n_y+1][n_x+1], P_I_lyL[n_y+1][n_x+1], U_I_lyL[n_y+1][n_x+1], V_I_lyL[n_y+1][n_x+1];
	double RHO_I_lyR[n_y+1][n_x+1], P_I_lyR[n_y+1][n_x+1], U_I_lyR[n_y+1][n_x+1], V_I_lyR[n_y+1][n_x+1];
	double RHO_F_lx[n_y+1][n_x+1], E_F_lx[n_y+1][n_x+1], U_F_lx[n_y+1][n_x+1], V_F_lx[n_y+1][n_x+1];
	double RHO_F_ly[n_y+1][n_x+1], E_F_ly[n_y+1][n_x+1], U_F_ly[n_y+1][n_x+1], V_F_ly[n_y+1][n_x+1];

  	double Z_g_MIDx[n_y+1][n_x+1], Z_g_MIDy[n_y+1][n_x+1];
  	double P_I_MIDx[n_y+1][n_x+1], U_I_MIDx[n_y+1][n_x+1], V_I_MIDx[n_y+1][n_x+1];
	double P_I_MIDy[n_y+1][n_x+1], U_I_MIDy[n_y+1][n_x+1], V_I_MIDy[n_y+1][n_x+1];

	double z_gL, z_gR, z_gxL, z_gxR, z_gyL, z_gyR;
	double rho_gL, rho_gR, p_gL ,p_gR ,u_gL ,u_gR ,v_gL ,v_gR;
	double rho_gxL, rho_gxR, p_gxL,p_gxR,u_gxL,u_gxR,v_gxL,v_gxR;
	double rho_gyL, rho_gyR, p_gyL,p_gyR,u_gyL,u_gyR,v_gyL,v_gyR;
	double rho_lL, rho_lR, p_lL ,p_lR ,u_lL ,u_lR ,v_lL ,v_lR;
	double rho_lxL, rho_lxR, p_lxL,p_lxR,u_lxL,u_lxR,v_lxL,v_lxR;
	double rho_lyL, rho_lyR, p_lyL,p_lyR,u_lyL,u_lyR,v_lyL,v_lyR;

	double z_gmid, rho_gmid, p_gmid, u_gmid, v_gmid;
	double z_lmid, rho_lmid, p_lmid, u_lmid, v_lmid;
	const double dx= config[10], dy= config[11];
	int iL,iR,jL,jR;
	double a_g, a_l, S_max, S_max_g, S_max_l;
	double wave_speed[2], dire[6], mid[6], star[6];

	double U[8];
	struct U_var U_L[n_y][n_x], U_R[n_y][n_x];
	double phi_sL[n_y][n_x], phi_sR[n_y][n_x];
	double S, S_tmp;
	for(int l = 0; l < (int)config[5] && stop_step != 1; ++l)
		{
			start_clock = clock();

			if (t_all >= plot_t)
				{
					for(int i = 0; i < n_y; ++i)
						for(j = 0; j < n_x; ++j)
							{
								FV->RHO[i*n_x+j]=ZRHO_gC[i][j]/Z_gC[i][j];
								FV->Z_a[i*n_x+j]=Z_gC[i][j];
								FV->PHI[i*n_x+j]=Z_gC[i][j];
							}
					file_write_TEC(*FV, *mv, problem, plot_t, dim);
					plot_t += delta_plot_t;
				}


			if (stop_step == 0)
				{									
					for(int i = 1; i < n_y; ++i)
						for(j = 1; j < n_x; ++j)
							{
								U[0] = ZRHO_gC[i][j];
								U[1] = RHO_U_gC[i][j];
								U[2] = RHO_V_gC[i][j];
								U[3] = E_gC[i][j]-0.5*RHO_V_gC[i][j]*RHO_V_gC[i][j]/ZRHO_gC[i][j];
								U[4] = ZRHO_lC[i][j];
								U[5] = RHO_U_lC[i][j];
								U[6] = RHO_V_lC[i][j];
								U[7] = E_lC[i][j]-0.5*RHO_V_lC[i][j]*RHO_V_lC[i][j]/ZRHO_lC[i][j];			
								phi_sL[i][j] = 0.5*(Z_gC[i-1][j-1]+Z_gC[i][j-1]);
								phi_sR[i][j] = 0.5*(Z_gC[i-1][j]+Z_gC[i][j]);					
								primitive_comp(U, &U_L[i][j], &U_R[i][j], phi_sL[i][j], phi_sR[i][j], 0.5, 0.5);
								RHO_gC[i][j] = U_L[i][j].lo_g;
								U_gC[i][j] = U_L[i][j].u_g;
								V_gC[i][j] = U_L[i][j].v_g;
								P_gC[i][j] = U_L[i][j].p_g;
								RHO_lC[i][j] = U_L[i][j].lo_s;
								U_lC[i][j] = U_L[i][j].u_s;
								V_lC[i][j] = U_L[i][j].v_s;
								P_lC[i][j] = U_L[i][j].p_s;
							}
				}
			else if (stop_step == 2)
				{									
					for(int i = 1; i < n_y; ++i)
						for(j = 1; j < n_x; ++j)
							{
								U[0] = ZRHO_gC[i][j];
								U[1] = RHO_V_gC[i][j];
								U[2] = RHO_U_gC[i][j];
								U[3] = E_gC[i][j]-0.5*RHO_U_gC[i][j]*RHO_U_gC[i][j]/ZRHO_gC[i][j];
								U[4] = ZRHO_lC[i][j];
								U[5] = RHO_V_lC[i][j];
								U[6] = RHO_U_lC[i][j];
								U[7] = E_lC[i][j]-0.5*RHO_U_lC[i][j]*RHO_U_lC[i][j]/ZRHO_lC[i][j];			
								phi_sL[i][j] = 0.5*(Z_gC[i-1][j-1]+Z_gC[i-1][j]);
								phi_sR[i][j] = 0.5*(Z_gC[i][j-1]+Z_gC[i][j]);					
								primitive_comp(U, &U_L[i][j], &U_R[i][j], phi_sL[i][j], phi_sR[i][j], 0.5, 0.5);
								RHO_gC[i][j] = U_L[i][j].lo_g;
								U_gC[i][j] = U_L[i][j].v_g;
								V_gC[i][j] = U_L[i][j].u_g;
								P_gC[i][j] = U_L[i][j].p_g;
								RHO_lC[i][j] = U_L[i][j].lo_s;
								U_lC[i][j] = U_L[i][j].v_s;
								V_lC[i][j] = U_L[i][j].u_s;
								P_lC[i][j] = U_L[i][j].p_s;
							}
				}			
			if (order > 1)
				{
					if (l==0)
						slope_limiter3_GRP(n_x, Z_gC, Z_I_gxL, Z_I_gxR, Z_I_gyL, Z_I_gyR, Z_gx, Z_gy, RHO_gC, P_gC, U_gC, V_gC, RHO_I_gxL, RHO_I_gxR, P_I_gxL, P_I_gxR, U_I_gxL, U_I_gxR, V_I_gxL, V_I_gxR, RHO_I_gyL, RHO_I_gyR, P_I_gyL, P_I_gyR, U_I_gyL, U_I_gyR, V_I_gyL, V_I_gyR, RHO_gx, P_gx, U_gx, V_gx, RHO_gy, P_gy, U_gy, V_gy, RHO_lC, P_lC, U_lC, V_lC, RHO_I_lxL, RHO_I_lxR, P_I_lxL, P_I_lxR, U_I_lxL, U_I_lxR, V_I_lxL, V_I_lxR, RHO_I_lyL, RHO_I_lyR, P_I_lyL, P_I_lyR, U_I_lyL, U_I_lyR, V_I_lyL, V_I_lyR, RHO_lx, P_lx, U_lx, V_lx, RHO_ly, P_ly, U_ly, V_ly);
					else
						slope_limiter_GRP(n_x, Z_gC, Z_I_gxL, Z_I_gxR, Z_I_gyL, Z_I_gyR, Z_gx, Z_gy, RHO_gC, P_gC, U_gC, V_gC, RHO_I_gxL, RHO_I_gxR, P_I_gxL, P_I_gxR, U_I_gxL, U_I_gxR, V_I_gxL, V_I_gxR, RHO_I_gyL, RHO_I_gyR, P_I_gyL, P_I_gyR, U_I_gyL, U_I_gyR, V_I_gyL, V_I_gyR, RHO_gx, P_gx, U_gx, V_gx, RHO_gy, P_gy, U_gy, V_gy, RHO_lC, P_lC, U_lC, V_lC, RHO_I_lxL, RHO_I_lxR, P_I_lxL, P_I_lxR, U_I_lxL, U_I_lxR, V_I_lxL, V_I_lxR, RHO_I_lyL, RHO_I_lyR, P_I_lyL, P_I_lyR, U_I_lyL, U_I_lyR, V_I_lyL, V_I_lyR, RHO_lx, P_lx, U_lx, V_lx, RHO_ly, P_ly, U_ly, V_ly);
				}

			if (stop_step == 0)
				{									
					tau = 1e15;
					for(int i = 1; i < n_y; ++i)
						for(j = 1; j < n_x; ++j)
							{
								a_g     = sqrt(gamma_g*P_gC[i][j]/RHO_gC[i][j]);
								a_l     = sqrt(gamma_l*P_lC[i][j]/RHO_lC[i][j]);
								S_max_g = fmax(fmax(fabs(U_gC[i][j]-a_g),fabs(U_gC[i][j]+a_g)),fmax(fabs(V_gC[i][j]-a_g),fabs(V_gC[i][j]+a_g)));
								S_max_l = fmax(fmax(fabs(U_lC[i][j]-a_l),fabs(U_lC[i][j]+a_l)),fmax(fabs(V_lC[i][j]-a_l),fabs(V_lC[i][j]+a_l)));
								S_max = fmax(S_max_g,S_max_l);
								tau   = fmin(tau,dx*config[7]/S_max);
							}
				}
			if((t_all + tau) > config[1])
				{
					printf("\nThe time is enough at step %d.\n",l);
					tau = config[1] - t_all;
					if (stop_step == 2)
						stop_t = 1;
				} // Time
			if (stop_step == 2)
				t_all += tau;

			if (stop_step == 0)
				{									
					for(int i = 1; i <  n_y; ++i)
						for(j = 1; j <= n_x; ++j)
							{
								if (j==1)
									{
										jL=1;
										jR=1;
									}
								else if (j==n_x)
									{
										jL=n_x-1;
										jR=n_x-1;
									}
								else
									{
										jL=j-1;
										jR=j;
									}
								z_gmid = U_R[i][jL].phi_s;
								z_lmid = 1.0-z_gmid;

								rho_gxL=RHO_gx[i][jL];
								rho_gxR=RHO_gx[i][jR];
								p_gxL=P_gx[i][jL];
								p_gxR=P_gx[i][jR];
								u_gxL=U_gx[i][jL];
								u_gxR=U_gx[i][jR];
								v_gxL=V_gx[i][jL];
								v_gxR=V_gx[i][jR];
								rho_gyL=RHO_gy[i][jL];
								rho_gyR=RHO_gy[i][jR];
								p_gyL=P_gy[i][jL];
								p_gyR=P_gy[i][jR];
								u_gyL=U_gy[i][jL];
								u_gyR=U_gy[i][jR];
								v_gyL=V_gy[i][jL];
								v_gyR=V_gy[i][jR];
								rho_gL =U_R[i][jL].lo_g+dx/2*rho_gxL;
								rho_gR =U_L[i][jR].lo_g-dx/2*rho_gxR;
								p_gL =U_R[i][jL].p_g+dx/2*p_gxL;
								p_gR =U_L[i][jR].p_g-dx/2*p_gxR;
								u_gL =U_R[i][jL].u_g+dx/2*u_gxL;
								u_gR =U_L[i][jR].u_g-dx/2*u_gxR;
								v_gL =U_R[i][jL].v_g+dx/2*v_gxL;
								v_gR =U_L[i][jR].v_g-dx/2*v_gxR;

								rho_lxL=RHO_lx[i][jL];
								rho_lxR=RHO_lx[i][jR];
								p_lxL=P_lx[i][jL];
								p_lxR=P_lx[i][jR];
								u_lxL=U_lx[i][jL];
								u_lxR=U_lx[i][jR];
								v_lxL=V_lx[i][jL];
								v_lxR=V_lx[i][jR];
								rho_lyL=RHO_ly[i][jL];
								rho_lyR=RHO_ly[i][jR];
								p_lyL=P_ly[i][jL];
								p_lyR=P_ly[i][jR];
								u_lyL=U_ly[i][jL];
								u_lyR=U_ly[i][jR];
								v_lyL=V_ly[i][jL];
								v_lyR=V_ly[i][jR];
								rho_lL =U_R[i][jL].lo_s+dx/2*rho_lxL;
								rho_lR =U_L[i][jR].lo_s-dx/2*rho_lxR;
								p_lL =U_R[i][jL].p_s+dx/2*p_lxL;
								p_lR =U_L[i][jR].p_s-dx/2*p_lxR;
								u_lL =U_R[i][jL].u_s+dx/2*u_lxL;
								u_lR =U_L[i][jR].u_s-dx/2*u_lxR;
								v_lL =U_R[i][jL].v_s+dx/2*v_lxL;
								v_lR =U_L[i][jR].v_s-dx/2*v_lxR;

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_gL, rho_gR, rho_gxL, rho_gxR, rho_gyL, rho_gyR, u_gL, u_gR, u_gxL, u_gxR, u_gyL, u_gyR, v_gL, v_gR, v_gxL, v_gxR, v_gyL, v_gyR, p_gL, p_gR, p_gxL, p_gxR, p_gyL, p_gyR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid[0] + 0.5*tau*dire[0];
								u_gmid   = mid[1] + 0.5*tau*dire[1];
								v_gmid   = mid[2] + 0.5*tau*dire[2];
								p_gmid   = mid[3] + 0.5*tau*dire[3];
								RHO_I_gxL[i][j] = mid[0] + tau*dire[0];
								U_I_gxL[i][j]   = mid[1] + tau*dire[1];
								V_I_gxL[i][j]   = mid[2] + tau*dire[2];
								P_I_gxL[i][j]   = mid[3] + tau*dire[3];

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_lL, rho_lR, rho_lxL, rho_lxR, rho_lyL, rho_lyR, u_lL, u_lR, u_lxL, u_lxR, u_lyL, u_lyR, v_lL, v_lR, v_lxL, v_lxR, v_lyL, v_lyR, p_lL, p_lR, p_lxL, p_lxR, p_lyL, p_lyR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_lmid = mid[0] + 0.5*tau*dire[0];
								u_lmid   = mid[1] + 0.5*tau*dire[1];
								v_lmid   = mid[2] + 0.5*tau*dire[2];
								p_lmid   = mid[3] + 0.5*tau*dire[3];
								RHO_I_lxL[i][j] = mid[0] + tau*dire[0];
								U_I_lxL[i][j]   = mid[1] + tau*dire[1];
								V_I_lxL[i][j]   = mid[2] + tau*dire[2];
								P_I_lxL[i][j]   = mid[3] + tau*dire[3];

								RHO_F_gx[i][j] = z_gmid*rho_gmid*u_gmid;
								U_F_gx[i][j]   = RHO_F_gx[i][j]*u_gmid + z_gmid*p_gmid;
								V_F_gx[i][j]   = RHO_F_gx[i][j]*v_gmid;
								E_F_gx[i][j]   = gamma_g/(gamma_g-1.0)*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								E_F_gx[i][j]   = RHO_F_gx[i][j]*E_F_gx[i][j];
								RHO_F_lx[i][j] = z_lmid*rho_lmid*u_lmid;
								U_F_lx[i][j]   = RHO_F_lx[i][j]*u_lmid + z_lmid*p_lmid;
								V_F_lx[i][j]   = RHO_F_lx[i][j]*v_lmid;
								E_F_lx[i][j]   = gamma_l/(gamma_l-1.0)*p_lmid/rho_lmid + 0.5*(u_lmid*u_lmid + v_lmid*v_lmid);
								E_F_lx[i][j]   = RHO_F_lx[i][j]*E_F_lx[i][j];

								Z_g_MIDx[i][j] = z_gmid;
								P_I_MIDx[i][j] = z_gmid*p_gmid + z_lmid*p_lmid;
								U_I_MIDx[i][j] = (z_gmid*rho_gmid*u_gmid + z_lmid*rho_lmid*u_lmid)/(z_gmid*rho_gmid + z_lmid*rho_lmid);
								V_I_MIDx[i][j] = (z_gmid*rho_gmid*v_gmid + z_lmid*rho_lmid*v_lmid)/(z_gmid*rho_gmid + z_lmid*rho_lmid);
								Z_I_gxR[i][j]   = Z_I_gxL[i][j];
								RHO_I_gxR[i][j] = RHO_I_gxL[i][j];
								U_I_gxR[i][j]   = U_I_gxL[i][j];
								V_I_gxR[i][j]   = V_I_gxL[i][j];
								P_I_gxR[i][j]   = P_I_gxL[i][j];
								RHO_I_lxR[i][j] = RHO_I_lxL[i][j];
								U_I_lxR[i][j]   = U_I_lxL[i][j];
								V_I_lxR[i][j]   = V_I_lxL[i][j];
								P_I_lxR[i][j]   = P_I_lxL[i][j];				
							}
					for(int i = 1; i < n_y; ++i)
						for(j = 1; j < n_x; ++j)
							{
								if(fabs(U_R[i][j].phi_s-U_L[i][j].phi_s)<eps)
									S = 0.5*(U_L[i][j].p_g+U_R[i][j].p_g)*(U_R[i][j].phi_s-U_L[i][j].phi_s);
								else
									S_tmp=U_R[i][j].phi_s*P_I_gxL[i][j+1]-U_L[i][j].phi_s*P_I_gxR[i][j];
									
								ZRHO_gC[i][j] -= tau*(RHO_F_gx[i][j+1]-RHO_F_gx[i][j])/dx;
								RHO_U_gC[i][j]-= tau*(U_F_gx[i][j+1]  -U_F_gx[i][j])  /dx-tau/dx*S;
								RHO_V_gC[i][j]-= tau*(V_F_gx[i][j+1]  -V_F_gx[i][j])  /dx;
								E_gC[i][j]    -= tau*(E_F_gx[i][j+1]  -E_F_gx[i][j])  /dx+tau/dx*S*U_L[i][j].u_s;
								RHO_U_lC[i][j]-= tau*(U_F_lx[i][j+1]  -U_F_lx[i][j])  /dx+tau/dx*S;
								RHO_V_lC[i][j]-= tau*(V_F_lx[i][j+1]  -V_F_lx[i][j])  /dx;
								E_lC[i][j]    -= tau*(E_F_lx[i][j+1]  -E_F_lx[i][j])  /dx+tau/dx*S*U_L[i][j].u_s;
							}

				}
			else if (stop_step == 2)
				{
					for(int i = 1; i <= n_y; ++i)
						for(j = 1; j <  n_x; ++j)
							{							
								if (i==1)
									{
										iL=1;
										iR=1;
									}
								else if (i==n_y)
									{
										iL=n_y-1;
										iR=n_y-1;
									}
								else
									{
										iL=i-1;
										iR=i;
									}
								z_gmid = U_R[iL][j].phi_s;
								z_lmid = 1.0-z_gmid;

								rho_gxL=RHO_gx[iL][j];
								rho_gxR=RHO_gx[iR][j];
								p_gxL=P_gx[iL][j];
								p_gxR=P_gx[iR][j];
								u_gxL=U_gx[iL][j];
								u_gxR=U_gx[iR][j];
								v_gxL=V_gx[iL][j];
								v_gxR=V_gx[iR][j];
								rho_gyL=RHO_gy[iL][j];
								rho_gyR=RHO_gy[iR][j];
								p_gyL=P_gy[iL][j];
								p_gyR=P_gy[iR][j];
								u_gyL=U_gy[iL][j];
								u_gyR=U_gy[iR][j];
								v_gyL=V_gy[iL][j];
								v_gyR=V_gy[iR][j];
								rho_gL =U_R[iL][j].lo_g+dy/2*rho_gyL;
								rho_gR =U_L[iR][j].lo_g-dy/2*rho_gyR;
								p_gL =U_R[iL][j].p_g+dy/2*p_gyL;
								p_gR =U_L[iR][j].p_g-dy/2*p_gyR;
								u_gL =U_R[iL][j].v_g+dy/2*u_gyL;
								u_gR =U_L[iR][j].v_g-dy/2*u_gyR;
								v_gL =U_R[iL][j].u_g+dy/2*v_gyL;
								v_gR =U_L[iR][j].u_g-dy/2*v_gyR;

								rho_lxL=RHO_lx[iL][j];
								rho_lxR=RHO_lx[iR][j];
								p_lxL=P_lx[iL][j];
								p_lxR=P_lx[iR][j];
								u_lxL=U_lx[iL][j];
								u_lxR=U_lx[iR][j];
								v_lxL=V_lx[iL][j];
								v_lxR=V_lx[iR][j];
								rho_lyL=RHO_ly[iL][j];
								rho_lyR=RHO_ly[iR][j];
								p_lyL=P_ly[iL][j];
								p_lyR=P_ly[iR][j];
								u_lyL=U_ly[iL][j];
								u_lyR=U_ly[iR][j];
								v_lyL=V_ly[iL][j];
								v_lyR=V_ly[iR][j];
								rho_lL =U_R[iL][j].lo_s+dy/2*rho_lyL;
								rho_lR =U_L[iR][j].lo_s-dy/2*rho_lyR;
								p_lL =U_R[iL][j].p_s+dy/2*p_lyL;
								p_lR =U_L[iR][j].p_s-dy/2*p_lyR;
								u_lL =U_R[iL][j].v_s+dy/2*u_lyL;
								u_lR =U_L[iR][j].v_s-dy/2*u_lyR;
								v_lL =U_R[iL][j].u_s+dy/2*v_lyL;
								v_lR =U_L[iR][j].u_s-dy/2*v_lyR;

								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_gL, rho_gR, rho_gyL, rho_gyR, -rho_gxL, -rho_gxR, v_gL, v_gR, v_gyL, v_gyR, -v_gxL, -v_gxR, -u_gL, -u_gR, -u_gyL, -u_gyR, u_gxL, u_gxR, p_gL, p_gR, p_gyL, p_gyR, -p_gxL, -p_gxR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_gmid = mid[0] + 0.5*tau*dire[0];
								u_gmid   =-mid[2] - 0.5*tau*dire[2];
								v_gmid   = mid[1] + 0.5*tau*dire[1];
								p_gmid   = mid[3] + 0.5*tau*dire[3];
								RHO_I_gyL[i][j] = mid[0] + tau*dire[0];
								U_I_gyL[i][j]   =-mid[2] - tau*dire[2];
								V_I_gyL[i][j]   = mid[1] + tau*dire[1];								
								P_I_gyL[i][j]   = mid[3] + tau*dire[3];
								linear_GRP_solver_Edir_Q1D(wave_speed, dire, mid, star, 0.0, 0.0, rho_lL, rho_lR, rho_lyL, rho_lyR, -rho_lxL, -rho_lxR, v_lL, v_lR, v_lyL, v_lyR, -v_lxL, -v_lxR, -u_lL, -u_lR, -u_lyL, -u_lyR, u_lxL, u_lxR, p_lL, p_lR, p_lyL, p_lyR, -p_lxL, -p_lxR, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, gamma_g, gamma_g, eps, eps);
								rho_lmid = mid[0] + 0.5*tau*dire[0];
								u_lmid   =-mid[2] - 0.5*tau*dire[2];
								v_lmid   = mid[1] + 0.5*tau*dire[1];
								p_lmid   = mid[3] + 0.5*tau*dire[3];
								RHO_I_lyL[i][j] = mid[0] + tau*dire[0];
								U_I_lyL[i][j]   =-mid[2] - tau*dire[2];
								V_I_lyL[i][j]   = mid[1] + tau*dire[1];								
								P_I_lyL[i][j]   = mid[3] + tau*dire[3];
								
								RHO_F_gy[i][j] = z_gmid*rho_gmid*v_gmid;
								U_F_gy[i][j]   = RHO_F_gy[i][j]*u_gmid;
								V_F_gy[i][j]   = RHO_F_gy[i][j]*v_gmid + z_gmid*p_gmid;
								E_F_gy[i][j]   = (gamma_g/(gamma_g-1.0))*p_gmid/rho_gmid + 0.5*(u_gmid*u_gmid + v_gmid*v_gmid);
								E_F_gy[i][j]   = RHO_F_gy[i][j]*E_F_gy[i][j];
								RHO_F_ly[i][j] = z_lmid*rho_lmid*v_lmid;
								U_F_ly[i][j]   = RHO_F_ly[i][j]*u_lmid;
								V_F_ly[i][j]   = RHO_F_ly[i][j]*v_lmid + z_lmid*p_lmid;
								E_F_ly[i][j]   = (gamma_l/(gamma_l-1.0))*p_lmid/rho_lmid + 0.5*(u_lmid*u_lmid + v_lmid*v_lmid);
								E_F_ly[i][j]   = RHO_F_ly[i][j]*E_F_ly[i][j];

								Z_g_MIDy[i][j] = z_gmid;
								P_I_MIDy[i][j] = z_gmid*p_gmid + z_lmid*p_lmid;
								U_I_MIDy[i][j] = (z_gmid*rho_gmid*u_gmid + z_lmid*rho_lmid*u_lmid)/(z_gmid*rho_gmid + z_lmid*rho_lmid);
								V_I_MIDy[i][j] = (z_gmid*rho_gmid*v_gmid + z_lmid*rho_lmid*v_lmid)/(z_gmid*rho_gmid + z_lmid*rho_lmid);

								Z_I_gyR[i][j]   = Z_I_gyL[i][j];
								RHO_I_gyR[i][j] = RHO_I_gyL[i][j];
								U_I_gyR[i][j]   = U_I_gyL[i][j];
								V_I_gyR[i][j]   = V_I_gyL[i][j];
								P_I_gyR[i][j]   = P_I_gyL[i][j];
								RHO_I_lyR[i][j] = RHO_I_lyL[i][j];
								U_I_lyR[i][j]   = U_I_lyL[i][j];
								V_I_lyR[i][j]   = V_I_lyL[i][j];
								P_I_lyR[i][j]   = P_I_lyL[i][j];
							}
					for(int i = 1; i < n_y; ++i)
						for(j = 1; j < n_x; ++j)
							{
								ZRHO_gC[i][j] -= tau*(RHO_F_gy[i+1][j]-RHO_F_gy[i][j])/dy;
								RHO_U_gC[i][j]-= tau*(U_F_gy[i+1][j]  -U_F_gy[i][j])  /dy;
								RHO_V_gC[i][j]-= tau*(V_F_gy[i+1][j]  -V_F_gy[i][j])  /dy;
								E_gC[i][j]    -= tau*(E_F_gy[i+1][j]  -E_F_gy[i][j])  /dy;
								RHO_U_lC[i][j]-= tau*(U_F_ly[i+1][j]  -U_F_ly[i][j])  /dy;
								RHO_V_lC[i][j]-= tau*(V_F_ly[i+1][j]  -V_F_ly[i][j])  /dy;
								E_lC[i][j]    -= tau*(E_F_ly[i+1][j]  -E_F_ly[i][j])  /dy;
							}
				}
						/*
						ZRHO_gC[i][j]  = Z_gC[i][j]*RHO_gC[i][j];
						RHO_U_gC[i][j] = U_gC[i][j]*ZRHO_gC[i][j];
						RHO_V_gC[i][j] = V_gC[i][j]*ZRHO_gC[i][j];
						E_gC[i][j] = Z_gC[i][j]*P_gC[i][j]/(gamma_g-1.0)+0.5*ZRHO_gC[i][j]*(U_gC[i][j]*U_gC[i][j]+V_gC[i][j]*V_gC[i][j]);
						ZRHO_lC[i][j]  = (1.0-Z_gC[i][j])*RHO_lC[i][j];
						RHO_U_lC[i][j] = U_lC[i][j]*ZRHO_lC[i][j];
						RHO_V_lC[i][j] = V_lC[i][j]*ZRHO_lC[i][j];
						E_lC[i][j] = (1.0-Z_gC[i][j])*P_lC[i][j]/(gamma_l-1.0)+0.5*ZRHO_lC[i][j]*(U_lC[i][j]*U_lC[i][j]+V_lC[i][j]*V_lC[i][j]);
						Z_gC[i][j]    -=tau*(0.5*(U_I_MIDx[i][j+1]+U_I_MIDx[i][j])*(Z_g_MIDx[i][j+1]-Z_g_MIDx[i][j])/dx+0.5*(V_I_MIDy[i+1][j]+V_I_MIDy[i][j])*(Z_g_MIDy[i+1][j]-Z_g_MIDy[i][j])/dy);
						ZRHO_gC[i][j] -=tau*((RHO_F_gx[i][j+1]-RHO_F_gx[i][j])/dx+(RHO_F_gy[i+1][j]-RHO_F_gy[i][j])/dy);
						RHO_U_gC[i][j]-=tau*((U_F_gx[i][j+1]  -U_F_gx[i][j])  /dx+(U_F_gy[i+1][j]  -U_F_gy[i][j])  /dy-0.5*(P_I_MIDx[i][j+1]+P_I_MIDx[i][j])*(Z_g_MIDx[i][j+1]-Z_g_MIDx[i][j])/dx);
						RHO_V_gC[i][j]-=tau*((V_F_gx[i][j+1]  -V_F_gx[i][j])  /dx+(V_F_gy[i+1][j]  -V_F_gy[i][j])  /dy-0.5*(P_I_MIDy[i+1][j]+P_I_MIDy[i][j])*(Z_g_MIDy[i+1][j]-Z_g_MIDy[i][j])/dy);
						E_gC[i][j]    -=tau*((E_F_gx[i][j+1]  -E_F_gx[i][j])  /dx+(E_F_gy[i+1][j]  -E_F_gy[i][j])  /dy-0.25*(P_I_MIDx[i][j+1]+P_I_MIDx[i][j])*(U_I_MIDx[i][j+1]+U_I_MIDx[i][j])*(Z_g_MIDx[i][j+1]-Z_g_MIDx[i][j])/dx-0.25*(P_I_MIDy[i+1][j]+P_I_MIDy[i][j])*(V_I_MIDy[i+1][j]+V_I_MIDy[i][j])*(Z_g_MIDy[i+1][j]-Z_g_MIDy[i][j])/dy);
						ZRHO_lC[i][j] -=tau*((RHO_F_lx[i][j+1]-RHO_F_lx[i][j])/dx+(RHO_F_ly[i+1][j]-RHO_F_ly[i][j])/dy);
						RHO_U_lC[i][j]-=tau*((U_F_lx[i][j+1]  -U_F_lx[i][j])  /dx+(U_F_ly[i+1][j]  -U_F_ly[i][j])  /dy+0.5*(P_I_MIDx[i][j+1]+P_I_MIDx[i][j])*(Z_g_MIDx[i][j+1]-Z_g_MIDx[i][j])/dx);
						RHO_V_lC[i][j]-=tau*((V_F_lx[i][j+1]  -V_F_lx[i][j])  /dx+(V_F_ly[i+1][j]  -V_F_ly[i][j])  /dy+0.5*(P_I_MIDy[i+1][j]+P_I_MIDy[i][j])*(Z_g_MIDy[i+1][j]-Z_g_MIDy[i][j])/dy);
						E_lC[i][j]    -=tau*((E_F_lx[i][j+1]  -E_F_lx[i][j])  /dx+(E_F_ly[i+1][j]  -E_F_ly[i][j])  /dy+0.25*(P_I_MIDx[i][j+1]+P_I_MIDx[i][j])*(U_I_MIDx[i][j+1]+U_I_MIDx[i][j])*(Z_g_MIDx[i][j+1]-Z_g_MIDx[i][j])/dx+0.25*(P_I_MIDy[i+1][j]+P_I_MIDy[i][j])*(V_I_MIDy[i+1][j]+V_I_MIDy[i][j])*(Z_g_MIDy[i+1][j]-Z_g_MIDy[i][j])/dy);
						*/						
			if (stop_step == 0)
				stop_step = 2;
			else if (stop_step == 2)
				stop_step = 0;			
			if (stop_t == 1)
				break;

			DispPro(t_all*100.0/config[1], l);
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;
		}

	for(int i = 0; i < n_y; ++i)
		for(j = 0; j < n_x; ++j)
			{
				FV->RHO[i*n_x+j]=ZRHO_gC[i][j]/Z_gC[i][j];
				FV->Z_a[i*n_x+j]=Z_gC[i][j];
				FV->PHI[i*n_x+j]=Z_gC[i][j];
			}

	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
