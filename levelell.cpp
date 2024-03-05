/*
 * levelell.cpp, is a implementation file for LevelEllipsoid class that 
 * generates the normal gravity field. 
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
 * Modified: 22.08.2010 A. Ustun    v2.1  added the second and third degree terms 
 *                                        of Taylor series of normal gravity 
 * 
 */

#include <stdio.h>
#include <math.h>
#include "levelell.h"

// constructor without parameter
LevelEllipsoid::LevelEllipsoid()
{
    /********** WGS84 **************/
    // a = 6378137.;
    //Rf = 298.257223563;
    //GM = 3986004.418e+8;
    // w = .7292115e-4;

    /********** GRS80 **************/
       a = 6378137.;
      J2 = 108263.0e-8;
      GM = 398600.5e+9;
       w = .7292115e-4;

    /********** GRS67 **************/
    // a = 6378160.;
    //J2 = 108270.0e-8;
    //GM = 398603.0e+9;
    // w = .72921151467e-4;

    /********** HAYFORD ************/
    // a = 6378388.;
    //Rf = 297.0;
    //ge = 9.78049;
    // w = .7292115e-4;

    double  comp = 0.0;
    double    q0 = 0.0;
    double   qu0 = 0.0;
    
    // computation of e1, e2 from physical parameters
    e1=sqrt(3*J2);
    do
    {
        comp=e1;
        e2=e1/sqrt(1-e1*e1);
        q0=0.5*((1+3/e2/e2)*atan(e2)-3/e2);
        qu0= 3*(1+1/e2/e2)*(1-atan(e2)/e2)-1;
        e1=sqrt(3*J2+4*w*w*a*a*a*e1*e1*e1/15/GM/2/q0);
    }while(fabs(e1-comp)>1.0e-13);
    
    
    // computation of geometric parameters
    Rf = 1/(1-sqrt(1-e1*e1));
    b  = a*(1-1/Rf);                      
    E  = sqrt(a*a-b*b);
    c  = a*a/b;            
    Q  = c*PI*(1-3*e2*e2/4+45*pow(e2,4)/64-175*pow(e2,6)/256+11025*pow(e2,8)/16384)/2;
    R1 = a*(1-1/Rf/3);
    R2 = c*(1-2*e2*e2/3+26*pow(e2,4)/45-100*pow(e2,6)/189+7034*pow(e2,8)/14175);
    R3 = pow(a*a*b,1/3.0);
    
    
    // computation of physical parameters
    m  = w*w*a*a*b/GM;
    U0 = GM*atan(e2)/E+w*w*a*a/3;
    ge = GM*(1-m-m*e2*qu0/6/q0)/a/b;
    gp = GM*(1+m*e2*qu0/3/q0)/a/a;
    kk = b*gp/a/ge-1.0;
    go = 1+e1*e1/6+kk/3+59*pow(e1,4)/360+5*e1*e1*kk/18
        +2371*pow(e1,6)/15120+259*pow(e1,4)*kk/1080+270229*pow(e1,8)/1814400+9623*pow(e1,6)*kk/45360;
    go*= ge;
    gf = (gp-ge)/ge;
    for(i=0; i<=10; i++)
        J[i]=pow(-1,i)*3*pow(e1,(i*2))*(1-i+5*i*J2/e1/e1)/(2*i+1)/(2*i+3);
}



// constructor with parameter
LevelEllipsoid::LevelEllipsoid(double pa, double pJ2,double pGM, double pw)
{
    a    = pa;
    J2   = pJ2;
    GM   = pGM;
    w    = pw;

    double  comp = 0.0;
    double    q0 = 0.0;
    double   qu0 = 0.0;

    if(J2<1.0)
    {
        // computation of e1, e2, Rf from physical parameters
        e1=sqrt(3*J2);
        do
        {
            comp=e1;
            e2=e1/sqrt(1-e1*e1);
            q0=0.5*((1+3/e2/e2)*atan(e2)-3/e2);
            qu0= 3*(1+1/e2/e2)*(1-atan(e2)/e2)-1;
            e1=sqrt(3*J2+4*w*w*a*a*a*e1*e1*e1/15/GM/2/q0);
        }while(fabs(e1-comp)>1.0e-12);
        Rf = 1/(1-sqrt(1-e1*e1));
        b  = a*(1-1/Rf);                      
        m  = w*w*a*a*b/GM;
    }
    else
    {
        Rf = J2;
        e1 = sqrt(2/Rf-1/Rf/Rf);
        e2 = e1/sqrt(1-e1*e1);
        q0 = 0.5*((1+3/e2/e2)*atan(e2)-3/e2);
        qu0= 3*(1+1/e2/e2)*(1-atan(e2)/e2)-1;
        b  = a*(1-1/Rf);                      
        m  = w*w*a*a*b/GM;
        J2 = e1*e1*(1-2*m*e2/15/q0)/3; 
    }
    
    
    // computation of geometric parameters
    E  = sqrt(a*a-b*b);
    c  = a*a/b;            
    Q  = c*PI*(1-3*e2*e2/4+45*pow(e2,4)/64-175*pow(e2,6)/256+11025*pow(e2,8)/16384)/2;
    R1 = a*(1-1/Rf/3);
    R2 = c*(1-2*e2*e2/3+26*pow(e2,4)/45-100*pow(e2,6)/189+7034*pow(e2,8)/14175);
    R3 = pow(a*a*b,1/3.0);
    
    
    // computation of physical parameters
    U0 = GM*atan(e2)/E+w*w*a*a/3;
    ge = GM*(1-m-m*e2*qu0/6/q0)/a/b;
    gp = GM*(1+m*e2*qu0/3/q0)/a/a;
    kk = b*gp/a/ge-1.0;
    go = 1+e1*e1/6+kk/3+59*pow(e1,4)/360+5*e1*e1*kk/18
        +2371*pow(e1,6)/15120+259*pow(e1,4)*kk/1080+270229*pow(e1,8)/1814400+9623*pow(e1,6)*kk/45360;
    go*= ge;
    gf = (gp-ge)/ge;
    J[0]=1.0;
    for(i=1; i<=10; i++)
        J[i]=pow(-1,i)*3*pow(e1,(i*2))*(1-i+5*i*J2/e1/e1)/(2*i+1)/(2*i+3);
}


// print values of system defining to stdout
void LevelEllipsoid::print()
{
    printf("NORMAL GRAVITY FIELD:\n");
    printf("=============================================================\n\n");
    printf("Defining parameters:\n\n");
    printf(" Semimajor axis                      a = %20.5f\n",a);
    printf(" Geocentric gravitational\n");
    printf(" constant                           GM = %20.11e\n", GM);
    printf(" Dynamic form factor                J2 = %20.11e\n", J2);
    printf(" Angular velocity                    w = %20.11e\n", w);
    printf("\n\n");
    printf("Derived geometric parameters:\n\n");
    printf(" semiminor axis                      b = %20.5f\n",b);
    printf(" polar radius of curv.               c = %20.5f\n",c);
    printf(" linear excentricity                 E = %20.5f\n",E);
    printf(" first excentricity                 e1 = %20.11e\n",e1);
    printf(" second excentricity                e2 = %20.11e\n",e2);
    printf(" flattaning                          f = %20.11e\n",1/Rf);
    printf(" reverse flattaning                1/f = %20.11f\n",Rf);
    printf(" meridian quadrant                (1/4)= %20.5f\n",Q);
    printf(" arithmetic mean radii of ellipsoid R1 = %20.5f\n",R1);
    printf(" radius of sphere of same surface   R2 = %20.5f\n",R2);
    printf(" radius of sphere of same volume    R3 = %20.5f\n",R3);
    printf("\n\n");
    printf("Derived physical parameters:\n\n");
    printf(" Normal potential at\n");
    printf(" ellipsoid surface                  U0 = %20.5f\n",U0);
    printf("                                     m = %20.11e\n",m);
    printf(" gravity at equator                    = %20.11f\n",ge);
    printf(" gravity at poles                      = %20.11f\n",gp);
    printf(" average gravity over the ellipsoid    = %20.11f\n",go);
    printf(" gravity flattenning                   = %20.11f\n",gf);
    printf("\n\n");
    printf("Zonal spherical harmonics:\n\n");
    for(i=1;i<=10;i++)
      printf("%4i%20.11e%20.11e\n",i*2,-J[i],J[i]/sqrt(4*i+1));
    printf("=============================================================\n\n\n");
}


double LevelEllipsoid::Vpot(double lat, double hgt)
// computes normal gravitational potential
{
    double slat=lat;
    double radi=hgt;
    ell2sph(&slat,&radi);

    double P[21];
    P[0]=1.0;
    P[1]=sin(slat);
    for(i=2;i<21;i++)
        P[i]=-(i-1)*P[i-2]/i+(2*i-1)*P[1]*P[i-1]/i;

    double U=0.0;
    for(i=1;i<11;i++)
        U+=pow(a/radi,2*i)*-J[i]*P[2*i];
    U=GM/radi*(1-U);

    return U;
}

double LevelEllipsoid::Cpot(double lat, double hgt)
// computes normal centrifugal potential
{
    double slat=lat;
    double radi=hgt;
    ell2sph(&slat,&radi);

    return 0.5*w*w*radi*radi*pow(cos(slat),2);
}


double LevelEllipsoid::Upot(double lat, double hgt)
// computes total gravity potential
{

    return Vpot(lat,hgt)+Cpot(lat,hgt);
}


double LevelEllipsoid::gonE(double lat)
// normal gravity on level ellipsoid
{
    double acos2 = a*cos(lat)*cos(lat);
    double bsin2 = b*sin(lat)*sin(lat);
    
    return (ge*acos2+gp*bsin2)/sqrt(a*acos2+b*bsin2);
}


double LevelEllipsoid::gath(double lat, double h)
// normal gravity at ellipsoidal height
{
    //return gonE(lat)*(1-2*(1+1/Rf+m-2*sin(lat)*sin(lat)/Rf)*h/a+3*h*h/a/a);
    return gonE(lat)+dgh1(lat)*h+dgh2(lat)*h*h/2+dgh3(lat)*h*h*h/6;
}

void LevelEllipsoid::ell2sph(double *lat,double *hgt)
// converts ellipsoidal coordinates to spherical ones
{
    double N=a/sqrt(1-e1*e1*pow(sin(*lat),2));
    double slat=atan(((1-e1*e1)*N+*hgt)*tan(*lat)/(N+*hgt));
    double radi=sqrt(pow((N+*hgt)*cos(*lat),2)+pow(((1-e1*e1)*N+*hgt)*sin(*lat),2));
    *lat=slat;
    *hgt=radi;
}


double LevelEllipsoid::dVdr(double lat, double hgt)
// computes the derivative of normal gravitational potential wtr radial cooordinate
{
    double slat=lat;
    double radi=hgt;
    ell2sph(&slat,&radi);

    double P[21];
    P[0]=1.0;
    P[1]=sin(slat);
    for(i=2;i<21;i++)
        P[i]=-(i-1)*P[i-2]/i+(2*i-1)*P[1]*P[i-1]/i;

    double dVdr=0.0;
    for(i=1;i<11;i++)
        dVdr+=(1+2*i)*pow(a/radi,2*i)*-J[i]*P[2*i];
    dVdr=-GM*(1-dVdr)/radi/radi;

    return dVdr;
}

double LevelEllipsoid::dVdN(double lat, double hgt)
// computes the derivative of normal gravitational potential wtr latitude
{
    double slat=lat;
    double radi=hgt;
    double d1P[21];
    double P[21];

    ell2sph(&slat,&radi);

      P[0]=1.0;
      P[1]=sin(slat);
    d1P[0]=0.0;
    d1P[1]=cos(slat);
    for(i=2;i<21;i++)
    {
          P[i]=-(i-1)*  P[i-2]/i+(2*i-1)*sin(slat)*P[i-1]/i;
        d1P[i]=-(i-1)*d1P[i-2]/i+(2*i-1)*cos(slat)*P[i-1]/i+(2*i-1)*sin(slat)*d1P[i-1]/i;
    //    printf("%20.11e%20.11e\n",P[i],d1P[i]);
    }

    double dVdN=0.0;
    for(i=1;i<11;i++)
        dVdN-=pow(a/radi,2*i)*-J[i]*d1P[2*i];
    dVdN*=GM/radi/radi;

    return dVdN;
}

double LevelEllipsoid::dCdr(double lat, double hgt)
// computes normal gravitational potential
{
    double slat=lat;
    double radi=hgt;
    ell2sph(&slat,&radi);

    return w*w*radi*cos(slat)*cos(slat);
}

double LevelEllipsoid::dCdN(double lat, double hgt)
// computes normal gravitational potential
{
    double slat=lat;
    double radi=hgt;
    ell2sph(&slat,&radi);

    return -w*w*radi*sin(slat)*cos(slat);
}

double LevelEllipsoid::gamm(double lat, double hgt)
// computes normal gravity at h using the norm of [grad(Vpot)+grad(Cpot)]
{
    double Vr=dVdr(lat,hgt); // radial   component of the derivative of gravitational pot
    double Vn=dVdN(lat,hgt); // latitude component of the derivative of gravitational pot
    double Cr=dCdr(lat,hgt); // radial   component of the derivative of centrifugal   pot
    double Cn=dCdN(lat,hgt); // latitude component of the derivative of centrifugal   pot
//    printf("%20.11e%20.11e%20.11e%20.11e\n",Vr,Vn,Cr,Cn);

    return sqrt(pow(Vr+Cr,2)+pow(Vn+Cn,2));
}

double LevelEllipsoid::Cgpu2Hnrm(double lat, double Cgpu)
// convert geopotential number to normal height
{
    double Hn=0.0;
    double g0=gonE(lat)/10;
    Hn=pow(Cgpu/a/g0,2);
    Hn+=1+Cgpu*(1+1/Rf+m-2*sin(lat)*sin(lat)/Rf)/a/g0;
    Hn*=Cgpu/g0;
    return Hn;
}

double LevelEllipsoid::W2Cgpu(double W, double K, double W0)
// convert geopotential number to normal height
{
    return (U0-W-(K-GM)/R1+W0-U0)/10.;
}

double LevelEllipsoid::Tzero(double K, double W0)
// convert geopotential number to normal height
{
    return (K-GM)/R1-(W0-U0);
}

double LevelEllipsoid::dgh1(double lat)
{
    double slat2=sin(lat)*sin(lat);
    return -gonE(lat)/a/(1-e1*e1)*sqrt(1-e1*e1*slat2)*(2-e1*e1-e1*e1*slat2)-2*w*w;
}
double LevelEllipsoid::dgh2(double lat)
{
    double slat2=sin(lat)*sin(lat);
    double carp=1-e1*e1*slat2;
    double a4=pow(a,4);
    return 6*GM/a4/carp/carp-30*GM*J2*(3*slat2-1)/carp/carp/carp/a4;
}
double LevelEllipsoid::dgh3(double lat)
{
    return -24*GM/pow(a,5);
}
