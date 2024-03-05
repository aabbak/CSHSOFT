/* comparison compares two data sets with the help of corrector surfaces. 
 * inputs : lat,long, data1, lat, lon, data2
 * outputs: min, max, mean, and rms of differences using 1, 4, 5, 7 parameter model
 *  Copyright (C) 2024 A. Abbak, KONYA, Türkiye
 *-----------------------------------------------------------------------------
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *-----------------------------------------------------------------------------
 * Author  : R. Alpay Abbak
 * Address : Konya Technical University 
 *           Geomatics Engineering Dept
 *           Kampus/Konya
 * E-Mail  : raabbak@ktun.edu.tr
 * Creation: 20.02.2010			v1.0 A. Abbak
 * Editing:  03.02.2024			v3.0 A. Abbak
 *
 * Compilation of the program on Linux: 
 * g++ comparison.cpp matris.cpp -o comparison
 *
 */
#include<stdio.h> 
#include<math.h>
#include"matris.h"
#define R 6371.0
#define PI 3.14159265359
#define mx 60000
#define rho 57.29577951308 
void HELP();
int main(int argc, char *argv[])
{/* Definitions ********************************************************************************************/
    FILE *stream1=NULL;					// first file pointer
    FILE *stream2=NULL;					// second file pointer 
    int 	   i=0;					// index
    int 	   j=0;					// index
    int           np=0;					// number of points
    int      baseN=0.0;					// Baseline number 
    double e2=0.006694380023;			        // the first eccentricity of WGS84
    double   sumv2=0.0;					// [vv]
    double    sumv=0.0;					// [v]
    double      m0=0.0;					// Root Mean Square Error 
    double     sph=0.0;					// Spherical Distance 
    double     ppm=0.0;					// parts per million 
    double     dif=0.0;					// differences
    double   total=0.0;					// total ppm 
    double    mean=0.0;					// mean ppm
    double     min=0.0;					// minimum ppm
    double     max=0.0;					// maximum ppm
    matris    phi1(mx);		 			// latitude of the first data
    matris    lam1(mx);					// longitude of the first data
    matris    phi2(mx);					// latitude of the second data
    matris    lam2(mx);					// longitude of the second data
    matris     geo(mx);					// Ngeometric
    matris     gra(mx);					// Ngravimetric
    matris       W(mx);					// coefficient 
/* File operations ***************************************************************************************/
    if(argc<3) 	HELP();    
    if((stream1=fopen(argv[1],"r")) == NULL)     // it consists of phi,lamda,Ngeometric
    {
		printf("The first file could not be opened!\n");
		exit(1);
    }
    if((stream2=fopen(argv[2],"r")) == NULL)     // it consists of phi,lamda,Ngravimetric
    {
		printf("The second file could not be opened!\n");
		exit(1);
    }
    while(!feof(stream1))
    {
	fscanf(stream1,"%lf%lf%lf\n",&phi1(i),&lam1(i),&geo(i));
	phi1(i)=phi1(i)/rho; lam1(i)=lam1(i)/rho;
	i++;
    }
    fclose(stream1);
    np=i;
    i=0;
    while(!feof(stream2))
    {
	fscanf(stream2,"%lf%lf%lf\n",&phi2(i),&lam2(i),&gra(i));
	if(fabs(phi1(i)-phi2(i)/rho)>1e-6 && fabs(lam1(i)-lam2(i)/rho)>1e-6){
	printf("\n\tThe datasets are not same order!!!\n"); exit(1);}
	i++;
    }
    fclose(stream2);
/* One Parameter model ***************************************************************************************/
   matris AAAA(np,1);
   matris NNNN(1,1);
   matris QQQQ(1,1);
   matris nnnn(1);
   matris xxxx(1);
   matris v(np);
   matris l(np);
   for(i=0;i<np;i++)
   {
	   AAAA(i,0)=1.0;					// design matrix
	   l(i,0)=(geo(i)-gra(i))*100.0;			// Observable vector ...
   }
   NNNN=transpoze(AAAA)*AAAA;
   QQQQ=invch(NNNN);
   nnnn=transpoze(AAAA)*l;
   xxxx=QQQQ*nnnn;
   v=AAAA*xxxx-l;
   sumv2=0.0;
   sumv=0.0;
   min=100.0;
   max=-100.0;
   mean=0.0;
   for(i=0;i<np;i++)
   {
	   if(min>l(i)) min=l(i);
	   if(max<l(i)) max=l(i);
	   sumv+=v(i);
	   sumv2+=v(i)*v(i);
   }
   m0=sqrt(sumv2/(np-1.0));
   printf("For 1 parameter model: Min=%6.2lf cm Max=%6.2lf cm Mean=%6.2lf cm m0=%6.2lf cm\n",min,max,xxxx(0),m0);
/* Four Parameter model ****************************************************************************************/
   matris A(np,4);
   matris N(4,4);
   matris Q(4,4);
   matris n(4);
   matris x(4);
  for(i=0;i<np;i++)
   {
	   A(i,0)=cos(phi1(i))*cos(lam1(i));		// Design Matrix filling ...
	   A(i,1)=cos(phi1(i))*sin(lam1(i));
	   A(i,2)=sin(phi1(i));
	   A(i,3)=1.0;
   }
   N=transpoze(A)*A;
   Q=invch(N);
   n=transpoze(A)*l;
   x=Q*n;
   v=A*x-l;
   sumv2=0.0;
   sumv=0.0;
   min=100.0;
   max=-100.0;
   mean=0.0;
   for(i=0;i<np;i++)
   {
	   if(min>v(i)) min=v(i);
	   if(max<v(i)) max=v(i);
	   sumv+=v(i);
	   sumv2+=v(i)*v(i);
   }
   m0=sqrt(sumv2/(np-4.0));
   printf("For 4 parameter model: Min=%6.2lf cm Max=%6.2lf cm Mean=%6.2lf cm m0=%6.2lf cm\n",min,max,1.0*sumv/np,m0);
/* Five Parameter model ****************************************************************************************/
   matris AA(np,5);
   matris NN(5,5);
   matris QQ(5,5);
   matris nn(5);
   matris xx(5);
   for(i=0;i<np;i++)
   {
	   AA(i,0)=cos(phi1(i))*cos(lam1(i));		// Design Matrix filling ...
	   AA(i,1)=cos(phi1(i))*sin(lam1(i));
	   AA(i,2)=sin(phi1(i));
	   AA(i,3)=1.0;
	   AA(i,4)=sin(phi1(i))*sin(phi1(i));
   }
   NN=transpoze(AA)*AA;
   QQ=invch(NN);
   nn=transpoze(AA)*l;
   xx=QQ*nn;
   v=AA*xx-l;
   sumv2=0.0;
   sumv=0.0;
   min=100.0;
   max=-100.0;
   mean=0.0;
   for(i=0;i<np;i++)
   {
	   if(min>v(i)) min=v(i);
	   if(max<v(i)) max=v(i);
	   sumv+=v(i);
	   sumv2+=v(i)*v(i);
   }
   m0=sqrt(sumv2/(np-5.0));
   printf("For 5 parameter model: Min=%6.2lf cm Max=%6.2lf cm Mean=%6.2lf cm m0=%6.2lf cm\n",min,max,1.0*sumv/np,m0);
/* Seven Parameter model ***************************************************************************************/
   matris AAA(np,7);
   matris NNN(7,7);
   matris QQQ(7,7);
   matris nnn(7);
   matris xxx(7);
   for(i=0;i<np;i++)
   {
	   W(i)=sqrt(1.0-e2*sin(phi1(i))*sin(phi1(i)));
	   AAA(i,0)=cos(phi1(i))*cos(lam1(i));	// Design Matrix filling ...
	   AAA(i,1)=cos(phi1(i))*sin(lam1(i));
	   AAA(i,2)=sin(phi1(i));
	   AAA(i,3)=cos(phi1(i))*sin(phi1(i))*cos(lam1(i))/W(i);
	   AAA(i,4)=cos(phi1(i))*sin(phi1(i))*sin(lam1(i))/W(i);
	   AAA(i,5)=sin(phi1(i))*sin(phi1(i))/W(i);
	   AAA(i,6)=1.0;
   }
   NNN=transpoze(AAA)*AAA;
   QQQ=invch(NNN);
   nnn=transpoze(AAA)*l;
   xxx=QQQ*nnn;
   v=AAA*xxx-l;
   sumv2=0.0;
   sumv=0.0;
   min=100.0;
   max=-100.0;
   mean=0.0;
   for(i=0;i<np;i++)
   {
	   if(min>v(i)) min=v(i);
	   if(max<v(i)) max=v(i);
	   sumv+=v(i);
	   sumv2+=v(i)*v(i);
   }
   m0=sqrt(sumv2/(np-7.0));
   printf("For 7 parameter model: Min=%6.2lf cm Max=%6.2lf cm Mean=%6.2lf cm m0=%6.2lf cm\n",min,max,1.0*sumv/np,m0);
/* Relative comparison before fit ******************************************************************************/
   min=100.0;
   max=-100.0;
   for(i=0;i<np;i++)
   {
   	for(j=i+1;j<np;j++)
	{
		dif=fabs(geo(i)-geo(j)-gra(i)+gra(j));
		sph=acos(sin(phi1(i))*sin(phi1(j))+cos(phi1(i))*cos(phi1(j))*cos(lam1(j)-lam1(i)))*R;
		ppm=dif*1000.0/sph;
		total+=ppm;
		if(min>ppm) min=ppm;
		if(max<ppm) max=ppm;
	}
   }
   baseN=np*(np-1)/2.0;
   printf("Relative comparison for %i baselines before fit: Min=%6.2lf ppm Max=%6.2lf ppm Mean=%6.2lf ppm\n",baseN,min,max,total/baseN*1.0);
/* Relative comparison after fit **********************************************************************************/
   min=1000.0;
   max=-1000.0;
   total=0.0;
   for(i=0;i<np;i++) // The data is being calibrated herein.
	gra(i)+=(A(i,0)*x(0)+A(i,1)*x(1)+A(i,2)*x(2)+A(i,3)*x(3))/100.0;
   for(i=0;i<np;i++)
   {
   	for(j=i+1;j<np;j++)
	{
		dif=fabs(geo(i)-geo(j)-gra(i)+gra(j));
		sph=acos(sin(phi1(i))*sin(phi1(j))+cos(phi1(i))*cos(phi1(j))*cos(lam1(j)-lam1(i)))*R;
		ppm=dif*1000.0/sph;
		total+=ppm;
		if(min>ppm) min=ppm;
		if(max<ppm) max=ppm;
	}
   }
   baseN=np*(np-1)/2.0;
   printf("Relative comparison for %i baselines  after fit: Min=%6.2lf ppm Max=%6.2lf ppm Mean=%6.2lf ppm\n",baseN,min,max,total/baseN);
/* Finish *********************************************************************************************************/
    return 0;
}
void HELP()
{
	fprintf(stderr,"\n     The comparison compares two datasets by corrector surfaces.\n\n");
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"     comparison  data1<file>  data2<file>\n\n");
	fprintf(stderr,"PARAMETERS:\n");
	fprintf(stderr,"     data1:      latitude, longitude, the first data.\n\n");
	fprintf(stderr,"                 It includes randomly data.\n\n");
	fprintf(stderr,"     data2:      latitude, longitude, the second data.\n\n");
	fprintf(stderr,"                 It includes randomly data.\n\n");
	fprintf(stderr,"AUTHOR:          Dr. R. Alpay ABBAK (raabbak@ktun.edu.tr)\n\n");
	exit(EXIT_FAILURE); 
}
