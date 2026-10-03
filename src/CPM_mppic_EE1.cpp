/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void CPM::mppic_EE1(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);
    
    press_lithostatic(p,a,pgc,s);
    
    pressure_gradient(p,a,pgc,s);
    
    ALOOP
    test(i,j,k) = dTz(i,j,k);
    
    // stress and cellSum update
    volfrac_update(p,pgc,s,P.X,P.Y,P.Z);
    stress_snider(p,pgc,s);
    stress_gradient(p,a,pgc,s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {


    // Save old position for bed-change calculation
         P.XRK1[n] = P.X[n];
         P.YRK1[n] = P.Y[n];
         P.ZRK1[n] = P.Z[n];

        // advec_mppic_step1(p, a, P, s, pturb,
        //            P.X, P.Y, P.Z, P.U, P.V, P.W,
         //           F, G, H, 0.5);

        // Velocity update 1
        // P.U[n] = (P.U[n] + p->dtsed*F)/(1.0 + p->dtsed*Dpx);
        // P.V[n] = (P.V[n] + p->dtsed*G)/(1.0 + p->dtsed*Dpy);
        // P.W[n] = (P.W[n] + p->dtsed*H)/(1.0 + p->dtsed*Dpz);
        
        // Position update
        // P.X[n] = P.X[n] + p->dtsed*P.U[n];
        // P.Y[n] = P.Y[n] + p->dtsed*P.V[n];
        // P.Z[n] = P.Z[n] + p->dtsed*P.W[n];


//------AddedByTina


advec_mppic_step1(p, a, P, s, pturb,
            P.X, P.Y, P.Z, P.U, P.V, P.W,
            F, G, H, 0.5);


// -------------------------------------------------------------
// Initiation of particle motion - Coulomb criterion
// Option A: criterion is applied only before first mobilization
// -------------------------------------------------------------

if(P.Mobile[n] < 0.5)
{
    // ---------------------------------------------------------
    // Hydrodynamic quantities for initiation of motion
    // sampled 0.01 m above the resting particle
    // ---------------------------------------------------------
    double z_ini = P.Z[n] + 0.01;

    // Fluid velocity at initiation sampling point
    double uf_ini =
        p->ccipol1c(a->u, P.X[n], P.Y[n], z_ini);

    double vf_ini =
        p->ccipol2c(a->v, P.X[n], P.Y[n], z_ini);

    double wf_ini =
        p->ccipol3c(a->w, P.X[n], P.Y[n], z_ini);

    // True particle-fluid relative velocity
    double Urel_ini = uf_ini - P.U[n];
    double Vrel_ini = vf_ini - P.V[n];
    double Wrel_ini = wf_ini - P.W[n];

    // Total 3-D slip magnitude for the drag coefficient
    double Urel_mag_ini =
        sqrt(Urel_ini*Urel_ini
           + Vrel_ini*Vrel_ini
           + Wrel_ini*Wrel_ini);

    // Solid volume fraction at initiation sampling point
    double Ts_ini =
        p->ccipol4a(Ts, P.X[n], P.Y[n], z_ini);

    // Drag coefficient/rate evaluated at initiation point
    double Dp_ini =
        drag_model(p, P.D[n], P.RO[n],
                   Urel_mag_ini, Ts_ini);

    // Horizontal drag acceleration
    double a_drag_x = Dp_ini * Urel_ini;
    double a_drag_y = Dp_ini * Vrel_ini;

    // Pressure gradients at initiation sampling point
    double dPx_ini =
        p->ccipol4a(dPx, P.X[n], P.Y[n], z_ini);

    double dPy_ini =
        p->ccipol4a(dPy, P.X[n], P.Y[n], z_ini);

    // Horizontal pressure-gradient acceleration
    double a_press_x = -dPx_ini / P.RO[n];
    double a_press_y = -dPy_ini / P.RO[n];

    // Total horizontal hydrodynamic driving acceleration
    double a_drive_x = a_drag_x + a_press_x;
    double a_drive_y = a_drag_y + a_press_y;

    double a_drive =
        sqrt(a_drive_x*a_drive_x
           + a_drive_y*a_drive_y);

    // Coulomb resistance based on submerged particle weight
    double a_res =
        p->Q30
        * (P.RO[n] - p->W1)/P.RO[n]
        * fabs(p->W22);

    // Initiation of motion
    if(a_drive > a_res)
        P.Mobile[n] = 1.0;
}


// -------------------------------------------------------------
// Particle motion
// -------------------------------------------------------------
if(P.Mobile[n] < 0.5)
{
    // Particle has not reached the initiation criterion
    P.U[n] = 0.0;
    P.V[n] = 0.0;
    P.W[n] = 0.0;

    // X, Y, Z remain unchanged
}
else
{
    // Normal MP-PIC velocity update
    P.U[n] =
        (P.U[n] + p->dtsed*F)
        /(1.0 + p->dtsed*Dpx);

    P.V[n] =
        (P.V[n] + p->dtsed*G)
        /(1.0 + p->dtsed*Dpy);

    P.W[n] =
        (P.W[n] + p->dtsed*H)
        /(1.0 + p->dtsed*Dpz);

    // Position update
    P.X[n] = P.X[n] + p->dtsed*P.U[n];
    P.Y[n] = P.Y[n] + p->dtsed*P.V[n];
    P.Z[n] = P.Z[n] + p->dtsed*P.W[n];
}





    }
    /*
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // advec 2
        advec_mppic_step2(p, a, P, s, pturb,
                    P.X, P.Y, P.Z, P.U, P.V, P.W,
                    F, G, H, 0.5);
                    
        //F=G=H=0.0;

        // Velocity update 2
        P.U[n] += 0.5*p->dtsed*F;
        P.V[n] += 0.5*p->dtsed*G;
        P.W[n] += 0.5*p->dtsed*H;

        // Position update
        P.X[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*p->dtsed*P.U[n];
        P.Y[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*p->dtsed*P.V[n];
        P.Z[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n];
    }*/

    boundcheck(p,2);
    bedchange_update(p,pgc,2);
    bedchange(p,a,pgc,s,2);

    // parallel transfer
    P.xchange(p, pgc,bedch,2);
}
