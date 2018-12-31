/*********************************************************************
/* File name: VIPLimiter.h
/* 
/* Abstract: The subroutin implements the VIP limiter for simulations of 2D flows on structured grids.
/*           Note, this is only a limiter to limit the velocity vector V=(u,v) for 2D flows. 
/* 
/*           The specific implementation is mainly based on section 2.1-2.3 of the following paper: 
/*           
/*           G.Luttwak, J.Falcovitz, Slope limiting for vectors, A novel vector limiting algorithm,
/*           Int. J. Numer. Meth. Fluids., 2011 (65)1365-1375.
/*
/* Date: Apr 11, 2018
/* Author: Jian Cheng @ IAPCM
***********************************************************************/
#ifndef __VIPLIMITER__H__
#define __VIPLIMITER__H__
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

/**
 * Subroutine of using VIP limiter for 2D velocity vector
 * 
 * \param[in] neigh_cell_num: number of neighbor cells 
 *                            (Note: only 2D face neighbor cells are used here, thus, 0=<neigh_cell_num<=4)
 *
 * \param[in] Vave: matrix for the average velocity vectors of neighbor cells
 *                  i.e. [u^0_ave,v^0_ave]
 *                       [u^1_ave,v^1_ave]
 *                       [u^2_ave,v^2_ave]
 *                       ...
 *
 * \param[in] V0: vector for the cell average velocity of the target cell
 *
 * \param[in|out] Vp: vector for the velocity of a given location in the target cell where VIP limiter needs to be applied
 *
 * \param[out] : the limiting coefficient lambda ( in [0, 1] ) for gradient vector
 *
 */
double useVIPLimiter(int neigh_cell_num, double Vave[][2], double* V0, double* Vp);

///////////////////////////////////////////////////
//some subroutines called by useViPLimiter...
///////////////////////////////////////////////////
double getTriArea(double x0, double y0, double x1, double y1, double xp, double yp);

void getPerpendFoot(double x0, double y0, double x1, double y1, double xc, double yc, double* pf);

bool obtuseAngle(double x0, double y0, double xa, double ya, double xb, double yb);

bool insideSegment(double x0, double y0, double x1, double y1, double xp, double yp);

double insectionPoint(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double* Vp);

bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp);

bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* Vp);

bool insideTriCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda);

bool insideQuadCH(std::vector<std::vector<double> >& CH, bool flag, double* V0, double* Vp, double& lambda);
///////////////////////////////////////////////////
#endif
