#ifndef _INITDTATA_H
#define _INITDTATA_H

#define Diaph3 (1.25) // Domain length
#define Tcell 400  // Number of computing cells in theta direction
#define Tcell_plot 100 // Output zoom
#define D_PLOT_T (0.001) // Output time interval
#define DATAOUT "../data_out" //data out folder


#define EPS (1e-9)
#define Alpha (0.) // GRP limiter parameter
#define LIMITER_CONF 2  /* LIMITER<0, add VIP limiter; LIMITER>0, only minmod limiter;
			   abs(LIMITER)=1, original minmod limiter; abs(LIMITER)=2, VIP-like minmod limiter */
#define Ncell 400 // Number of computing cells in r direction
#define Diaph1 (1.1)
#define Diaph2 (1.1)
#define Domlen (2.) // Domain length
#define Timeout (0.1) // Output time

#define GAMMAL (1.6667)
#define GAMMAR (1.4) // Ratio of special heats Gamma=1.4 or 3.0
#define DL0 (9995.87681795254)
#define DM0 (2500.)
#define DR0 (0.1)
#define UL0 (1.73184512961398)
#define UR0 (0.)
#define PL0 (10000.)
#define PR0 (1.)


#endif
