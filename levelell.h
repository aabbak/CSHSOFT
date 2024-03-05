/*
 * levelell.h, v1.0, is a header file for LevelEllipsoid class that generates
 * the normal gravity field. 
 *
 * Copyright (C) 2006 A. Ustun
 * 
 *-----------------------------------------------------------------------------
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *-----------------------------------------------------------------------------
 *
 * Author  : Aydin Ustun
 * Address : Selcuk Universitesi
 *           Jeodezi ve Fotogrametri Muh.
 *           Kampus/Konya
 * E-Mail  : austun@selcuk.edu.tr
 * 
 * Created : 13.07.2005             v1.0   
 * Modified: 10.05.2006 A. Ustun    v1.1  Regularization of the source code
 *                                        for separate compilation
 * Modified: 16.03.2009 A. Ustun    v2.0  centrifugal+gravity potential member
 *                                        functions added
 *
 */

#ifndef LEVELELL_H
#define LEVELELL_H

#define PI (2*asin(1.0))

#include <stdio.h>

class LevelEllipsoid
{
  public:
    // constructors
    LevelEllipsoid();
    LevelEllipsoid(double, double, double, double);

    double flat(){return 1/Rf;}
    double gave(){return go;}
    double Rave(){return R2;}
    
    double Vpot(double, double);
    double Cpot(double, double);
    double Upot(double, double);
    double gonE(double);
    double gath(double, double);
    double sphcoefs(int i){return J[i];}
    double dVdr(double, double);
    double dVdN(double, double);
    double dCdr(double, double);
    double dCdN(double, double);
    double gamm(double, double);
    double dgh1(double); // the first derivative of normal gravity 
    double dgh2(double); // the second derivative of normal gravity 
    double dgh3(double); // the third derivative of normal gravity 
    // convert ellipsoidal latitude and height to spherical ones
    void ell2sph(double *,double *);
    // convert geopotential number to normal height
    double Cgpu2Hnrm(double, double);
    // convert gravity potential to geopotential number
    double W2Cgpu(double,double,double);
    // compute zero degree term in geopotential unit
    double Tzero(double,double);
    void print();                                  // prints values of system defining to stdout

  private:
    int i;
    
    double     a;         // semimajor axis
    double    J2;         // dynamical form factor
    double    GM;         // GM for normal gravity field
    double     w;         // angular velocity
    
    
    /******************* Derived geometric parameters *******************/
    double    e1;         // first  eccentricity of level ellipsoid
    double    e2;         // second eccentricity
    double    Rf;         // reverse flattening
    double     b;         // semiminor axis
    double     E;         // linear eccentricity
    double     c;         // radii of polar curvature
    double     Q;         // meridian quadrant
    double    R1;         // arithmetic mean radii of ellipsoid
    double    R2;         // radius of sphere of the same surface
    double    R3;         // radius of sphere of the same volume
    
    
    /******************* Derived physical parameters ********************/
    double     m;         // w**2a**2b/GM
    double    U0;         // normal potential at level ellipsoid
    double    ge;         // normal gravity at equator
    double    gp;         // normal gravity at poles
    double    kk;         // constant for gravity calculations
    double    go;         // average value of gravity over the ellipsoid
    double    gf;         // gravity flattening
    double J[11];         // zonal spherical harmonic coefficients  
    
};

#endif

