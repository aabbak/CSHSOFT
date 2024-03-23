/* The CSHSOFT computes a gravimetric geoid model by using classical Stokes-Helmert method.
 *
 * This program is distributed under the terms of a special License published by the authors.
 * Hopefully, this program will be useful, but WITHOUT ANY WARRANTY. See the License of the 
 * software for more details. 
 *
 * Reference paper for the program:
 * Abbak, R. A., Goyal, R. and Ustun, A. (2024), A user-friendly software package for computing 
 * geoid model by classical Stokes-Helmert method, Computers & Geosciences.
 *
 * Compilation of the program on Linux: 
 * g++ CSHSOFT.cpp matris.cpp levelell.cpp -o CSHSOFT
 *
 * Execution of the program with our sample data and default parameters: 
 * ./CSHSOFT -GXGM2019.gfc -Aanomaly.xyz -Eelevation.xyz -Ttc.xyz -M630  
 * 
 * Created : 23.03.2023	R. A. Abbak	v1.0  						*/
#include <unistd.h>				// Standard option library 
#include "levelell.h"				// User-defined level ellipsoid library
#include "matris.h"				// User-defined matrix library
#define G 6.6742e-11				// Newtonian Earth's attraction constant
#define GM 0.3986005e+15			// Earth mass*gravity constant
#define omega 7292115.0e-11			// Angular velocity of the Earth's rotation
#define	a 6378137.0				// Semimajor axes of GRS80 ellipsoid
#define J2 108263.0e-8				// Dynamical form factor of the Earth (GRS80)
#define R 6371000.79				// Earth's mean radius
#define pi 3.141592653589793			// constant pi
#define rho 57.2957795130823			// 180degree/pi
#define mx 800					// Maximum dimension of the data area 
#define OPTIONS "G:A:E:T:M:L:P:R:I:SH"		// Parameters and Options
void HELP();					// Help function for proper usage
void LEGENDRE(matris &P,matris &P1,int N,double x);// Fully normalized Legendre function
int main(int argc, char *argv[]){
	FILE 	 *model=NULL;			// GGM file
	FILE   *anomaly=NULL;			// anomaly data file
	FILE *elevation=NULL;			// elevation data file
	FILE   *terrain=NULL;			// terrain corr. model file
	char 	      option;			// temprorary option
	char       line[256];  			// read a line from the GGM file
    	const char *slash="/";			// discriminant symbol for options
	int 	         s=0;			// shows whether all segments are printed
    	int              i=0;			// array element
	int      imax=300000;			// maximum array element for computation
    	int              j=0;			// array element
        int              n=0;			// array element of degree
    	int              m=0;			// array element of order
	int           mmax=0;			// elements of max expansion of GGM used
    	int             in=0;			// array element of computation point
    	int             jn=0;			// array element of computation point
	int            M=630;			// maximum expansion of GGM used
	int            L=145;			// maximum expansion of Stokes series
	double psi0=0.95/rho;			// max capsize in radian
	double 	MinPhi=45.01;			// minimum latitude of target area
	double 	MaxPhi=47.00;			// maximum latitude of target area
	double 	 MinLam=1.51;			// minimum longitude of target area
	double 	 MaxLam=4.50;			// maximum longitude of target area
	double 	 PhiInt=0.02;			// grid size in latitude direction
	double 	 LamInt=0.02;			// grid size in longitude direction
	double MinDatPhi=0.0;			// minimum latitude of data area
	double MinDatLam=0.0;			// minimum longitude of data area
	double        HC=0.0; 			// Cnm coefficient/height of computation point
    	double        HS=0.0; 			// Snm coefficient
    	double        Hc=0.0; 			// variance of Cnm coefficient
    	double        Hs=0.0; 			// variance of Snm coefficient
	double       phi=0.0;			// latitude of computational point [deg] 
    	double       lam=0.0;			// longitude of computational point [deg]
	double         t=0.0;			// cos(psi0)
	double        tt=0.0;			// sin(psi0/2)
	double    total1=0.0;			// total value
    	double    total2=0.0;			// total value
	matris       C(imax);			// Cnm coefficients
	matris       S(imax);			// Snm coefficients
	matris 	    lati(mx);			// latitude of running point in radian
	matris      loni(mx);			// longitude of running point in radian
	matris     Dg(mx,mx);			// terrestrial gravity anomalies
	matris      H(mx,mx);			// topographic heights
	matris    DgH(mx,mx);			// terrain corection on gravity anomaly
	matris 	    PN(imax);			// UNnormalized Legendre Function
	matris 	     P(imax);			// Normalized Legendre Function
	matris 	    P1(imax);			// The first derivation of Legendre Function
	matris     gamma(mx);			// normal gravity on ellipsoid
	matris   NGGM(mx,mx);			// geoid height derived from GGM
	matris  DgGGM(mx,mx);			// gravity anomaly derived from GGM
	matris  DgEll(mx,mx);			// Ellipsoidal Effect
	matris  DgDAE(mx,mx);			// Direct Atmospheric Effect
	matris DgSITE(mx,mx);			// Secondary Indirect Topographic Effect
	matris  Dgred(mx,mx);			// reduced gravity anomaly
/******** O P T I O N S   A N A L Y S I S **************************************************************************/
	if(argc<5) HELP();
	while((option=getopt(argc,argv,OPTIONS))!=-1)
        switch(option){
		case 'G': model=fopen(optarg,"r"); break; 	// GGM file opens	
		case 'A': anomaly=fopen(optarg,"r"); break;	// Anomaly file opens
		case 'E': elevation=fopen(optarg,"r"); break;	// elevation file opens
		case 'T': terrain=fopen(optarg,"r"); break;	// terrain cor. file opens
		case 'M': M=atoi(optarg); break;		// maximum expansion GGM used
		case 'L': L=atoi(optarg); break;		// maximum expansion of spherical harmonic
		case 'P': psi0=atof(optarg)/rho;break;		// capsize
		case 'R': MinPhi=atof(strtok(optarg,slash));    // Limits of target area 
			  MaxPhi=atof(strtok(NULL,slash));
			  MinLam=atof(strtok(NULL,slash));
			  MaxLam=atof(strtok(NULL,slash)); break;
		case 'I': PhiInt=atof(strtok(optarg,slash));    // resolution of geoid model
                	  LamInt=atof(strtok(NULL,slash)); break;
		case 'S': s=1; break;				// all segments will be printed
		case 'H': HELP();
		default : HELP();
        }
/******** R E A D I N G   F I L E S *********************************************************************************/
	if(model==NULL){
		printf("\nGGM file cannot be opened!!!\n\n"); exit(EXIT_FAILURE);
    	}
	if(anomaly==NULL){
		printf("\nAnomaly file cannot be opened!!!\n\n"); exit(EXIT_FAILURE);
    	}
	if(elevation==NULL){
		printf("\nElevation file cannot be opened!!!\n\n"); exit(EXIT_FAILURE);
    	}
	if(terrain==NULL){
		printf("\nTerrain corr. file cannot be opened!!!\n\n"); exit(EXIT_FAILURE);
    	}
	while(fgets(line,256,model)!=NULL){
        	if(sscanf(line,"%i%i%lf%lf%lf%lf",&n,&m,&HC,&HS,&Hc,&Hs)==6){
			i=n*(n+1)/2+m;		C(i)=HC;	S(i)=HS;	
        	}
    	}
    	fclose(model);
	mmax=(M+1)*(M+2)/2;
	MinDatPhi=MinPhi-5.0;
	MinDatLam=MinLam-5.0;
	while(!feof(anomaly)){
		fscanf(anomaly,"%lf%lf",&phi,&lam);
		i=round((phi-MinDatPhi)/PhiInt); 	
		j=round((lam-MinDatLam)/LamInt);
		lati(i)=phi/rho; 	
		loni(j)=lam/rho;
		fscanf(anomaly,"%lf\n",&Dg(i,j));
    	}
    	fclose(anomaly);
	while(!feof(elevation)){
		fscanf(elevation,"%lf%lf",&phi,&lam);
		i=round((phi-MinDatPhi)/PhiInt); 	
		j=round((lam-MinDatLam)/LamInt);
		fscanf(elevation,"%lf\n",&H(i,j));
    	}
    	fclose(elevation);
	while(!feof(terrain)){
		fscanf(terrain,"%lf%lf",&phi,&lam);
		i=round((phi-MinDatPhi)/PhiInt); 	
		j=round((lam-MinDatLam)/LamInt);
		fscanf(terrain,"%lf\n",&DgH(i,j));
    	}
    	fclose(terrain);
/******** M O D I F I C A T I O N  C O E F F I C I E N T S ***********************************************************/
	double   Q0=0.0;					// Molodenski coefficient
	double  QL0=0.0;					// Molodenski coefficient
	matris E(L+2,0);					// Paul's coefficient
	PN(0)=1.0;
	PN(1)=t=cos(psi0);
	for(n=2;n<L+2;n++)					// Unnormalized Legendre function
		PN(n)=((2.0*n-1.0)/n)*t*PN(n-1)-((n-1.0)/n)*PN(n-2);  // ok
	for(n=2;n<=L;n++)
	{
 		E(n,0)=(PN(n+1)-PN(n-1))/(2.0*n+1.0);
		total2+=(2.0*n+1.0)/(n-1.0)*E(n,0); 		// sn was included
	}
	t=sin(psi0/2.0);
	Q0=-4.0*t+5.0*t*t+6.0*t*t*t-7.0*t*t*t*t+(6.0*t*t-6.0*t*t*t*t)*log((t)*(1.0+t));
	QL0=Q0-total2;
	total2=0.0;
/******** R E M O V E   S E C T I O N *********************************************************************************/
	double constant3=-2.0*pi*G*2670.0/R*1.0e+5; 		// coefficient for the SITE 
	double constant4=-1.0e+5*0.006694380022903*GM/R/R; 	// coefficient for elliposidal correction 
	double     phii=0.0;					// latitude of computational point [rad]
	double     slat=0.0;					// spherical latitude of computational point [rad]
    	double     radi=0.0;					// spherical radius
	double   coslon=0.0;					// cos(longitude)
    	double   sinlon=0.0;					// sin(longitude)
    	double      Rr1=0.0;					// R/r
    	double      Rrn=0.0;					// R/r^n
    	double       Yn=0.0;					// Ymn
    	double       Y1=0.0;					// Ymn
    	double     dpot=0.0;  					// disturbing potential [m**2/s**2]
    	double    dpot1=0.0;  					// disturbing potential [m**2/s**2]
    	double     gdst=0.0;  					// gravity disturbance [m/s**2]
	matris     Dgn(M+2); 					// Dg[n]
	matris cosmlon(M+2);					// cos(mLongitude)
  	matris sinmlon(M+2);					// sin(mLongitude)	
    	if(C(0)==0.0)       C(0)=1.0;
	LevelEllipsoid levell(a,J2,GM,omega);			// Create normal gravity field (GRS80) 
	for(i=0;i<11;i++)					// Compute harmonic coefficients of disturbing potential 
		C(2*i*i+i)-=levell.sphcoefs(i)/sqrt(4*i+1);
	cosmlon(0)=1.0;						// cos(longitude)
	for(phi=MinPhi-1.0;phi<=MaxPhi+1.0;phi=phi+PhiInt)
	{
		phii=phi/rho;
		slat=phii;
		radi=0.0;
		levell.ell2sph(&slat,&radi);
		LEGENDRE(P,P1,M,sin(slat));
		in=round((phi-MinDatPhi)/PhiInt);
		gamma(in)=9.7803267715*(1.0+0.001931851353*pow(sin(phii),2))/sqrt(1.0-0.006694380023*pow(sin(phii),2));		
		for(lam=MinLam-1.5;lam<=MaxLam+1.5;lam=lam+LamInt)
		{
			jn=round((lam-MinDatLam)/LamInt);
			cosmlon(1)=coslon=cos(lam/rho);
			sinmlon(1)=sinlon=sin(lam/rho);
			m=1;
			while(++m<=M)
			{
				cosmlon(m)=2.0*coslon*cosmlon(m-1)-cosmlon(m-2);
				sinmlon(m)=2.0*coslon*sinmlon(m-1)-sinmlon(m-2);
			}
			Rr1=a/radi;
			Rrn=1.0;
			n=2; m=0; i=3;
			while(i<mmax)
			{
				Yn +=P(i)*(C(i)*cosmlon(m)+S(i)*sinmlon(m));
				Y1 +=P1(i)*(C(i)*cosmlon(m)+S(i)*sinmlon(m));
				if(m==n)
				{
					Rrn*=Rr1;
					dpot+=Yn*Rrn;
					dpot1+=Y1*Rrn;
        				gdst+=Yn*Rrn*(n+1.0);
	       				Dgn(n)=GM/a/radi*(gdst-2.0*dpot);
					Yn=0.0;	Y1=0.0; n++; m=0;
				}
				else    m++;
				i++;
			}
	      		DgGGM(in,jn)=1.0e+5*Dgn(n-1); 		// ok
 			NGGM(in,jn)=dpot*GM/R/gamma(in);	// NGGM(in,jn)=R*dpot;
			DgEll(in,jn)=constant4*(sin(slat)*cos(slat)*dpot1-(3.0*sin(slat)*sin(slat)-2.0)*dpot);
			HC=H(in,jn);
			DgDAE(in,jn)=0.871-1.0298e-4*HC+5.3105e-9*HC*HC-2.1642e-13*HC*HC*HC+9.5246e-18*pow(HC,4)-2.2411e-22*pow(HC,5);// ok
			DgSITE(in,jn)=constant3*HC*HC;				 						// ok
			Dgred(in,jn)=Dg(in,jn)-DgGGM(in,jn)+DgH(in,jn)+DgSITE(in,jn)+DgDAE(in,jn)+DgEll(in,jn); 	     	// ok
//		 	printf("%.4f %.4f %9.4f %9.4f %9.4f %9.4f\n",phi,lam,DgGGM(in,jn),DgDAE(in,jn),DgSITE(in,jn),DgEll(in,jn));
			dpot=dpot1=gdst=0.0;
		}
	}
/******** C O M P U T E   S E C T I O N ***********************************************************************************************/
	int framePhi=round(psi0*rho/PhiInt);				// vertical limit of compartment
    	int frameLam=round(acos((cos(psi0)-pow(sin(MaxPhi/rho),2))/(pow(cos(MaxPhi/rho),2)))*rho/LamInt);//horizontal limit of compartment	
	int              x=0;						// array element of frame
    	int              y=0;						// array element of frame
	double constant0=R/4.0/pi*PhiInt*LamInt/rho/rho;		// grid size of the block 
	double constant1=-G*2670.0*R*R/6.0*PhiInt*LamInt/rho/rho;	// coefficient for the first part of NH 
	double constant2=-pi*G*2670.0; 					// coefficient for the second part of NH 
	double        St=0.0;						// Stokes kernel
	double        N2=0.0;						// N2 part of Stokes' integration
	matris    psi(mx,mx);						// computational psi between P and Q
	matris     l0(mx,mx);						// spherical distance between P and Q
	matris     SM(mx,mx);						// modified Stokes' function
	matris stokes(mx,mx);						// Stokes' integration part of the geoid
	matris     NH(mx,mx);						// topographic correction (PITE)	
	matris      N(mx,mx);						// final geoid height
	for(phi=MinPhi;phi<=MaxPhi;phi=phi+PhiInt)
	{
		lam=MinLam;
		in=round((phi-MinDatPhi)/PhiInt);
		jn=round((lam-MinDatLam)/LamInt);
		for(i=in-framePhi;i<=in+framePhi;i++)
		{
			x=i-in+framePhi; 				// position on the compartment
			for(j=jn;j<=jn+frameLam;j++)
			{
				if(i==in && j==jn) continue;
				y=abs(j-jn); 				// position on the compartment
				psi(x,y)=acos(sin(lati(in))*sin(lati(i))+cos(lati(in))*cos(lati(i))*cos(loni(j)-loni(jn)));
				if(psi(x,y)<psi0)
				{	
					tt=sin(psi(x,y)/2.0);
					l0(x,y)=pow(2.0*R*tt,3); 	// l0^3 computed
					PN(1)=t=cos(psi(x,y));	
					St=1.0/tt+1.0-5.0*t-6.0*tt-3.0*t*log(tt+tt*tt); // Stokes kernel
					for(n=2;n<=L;n++)
					{
						PN(n)=((2.0*n-1.0)/n)*t*PN(n-1)-((n-1.0)/n)*PN(n-2);
						total1+=(2.0*n+1.0)/(n-1.0)*PN(n);
					}
					SM(x,y)=(St-total1)*cos(lati(i));
					total1=0.0;
				}
			}
		}
		for(lam=MinLam;lam<=MaxLam;lam=lam+LamInt)
		{
			jn=round((lam-MinDatLam)/LamInt);
			for(i=in-framePhi;i<=in+framePhi;i++)
			{
				x=i-in+framePhi; 			// position on the compartment
				for(j=jn-frameLam;j<=jn+frameLam;j++)
				{
					if(i==in && j==jn) continue;
					y=abs(j-jn); 			// position on the compartment
					if(psi(x,y)<psi0)		// both computations are done in same capsize!!
					{
 						total1+=SM(x,y)*(Dgred(i,j)-Dgred(in,jn))*cos(lati(i)); // modified Stokes kernel
						total2+=(pow(H(i,j),3)-pow(H(in,jn),3))/l0(x,y)*cos(lati(i)); // second term of topographic effect
					}
				}
			}
			stokes(in,jn)=constant0*total1/gamma(in)/1.0e+5;
			N2=-R/2.0/gamma(in)*Dgred(in,jn)*QL0/1.0e+5;
			NH(in,jn)=(constant1*total2+constant2*H(in,jn)*H(in,jn))/gamma(in); 
/******** R E S T O R E   S E C T I O N **************************************************************************************/
			N(in,jn)=NGGM(in,jn)+stokes(in,jn)+N2+NH(in,jn);
		 	if(s) printf("%.4f %.4f %9.4f %9.4f %9.4f %9.4f\n",phi,lam,NGGM(in,jn),stokes(in,jn)+N2,NH(in,jn),N(in,jn));
			else printf("%.4f %.4f %9.4f\n",phi,lam,N(in,jn));
			total1=total2=0.0;
		}
	}
	return 0;
/******** F I N I S H *******************************************************************************************************/
}
void LEGENDRE(matris &P, matris &P1, int N, double x){
	int i=0;
    	int j=0;
    	int k=0;
    	int n=0;
    	int m=0;
    	int imax=(N+1)*(N+2)/2;
    	long double sinx=x; 
    	long double cosx=sqrtl(1-x*x);
    	long double f1=.0;
    	long double f2=.0;
    	long double f3=.0;
    	long double f4=.0;
    	long double f5=.0;
    	P(0)=1.0;
    	P(1)=sqrt(3.0)*sinx;
    	P(2)=sqrt(3.0)*cosx;
    	P1(0)=0.0;
    	P1(1)= sqrtl(3.0)*cosx;
    	P1(2)=-sqrtl(3.0)*sinx;
    	for(n=2;n<=N;n++){
    		i=n*(n+1)/2+n; 						// index for Pn,n
        	j=i-n-1;       						// index for Pn-1,n-1
        	f1=sqrtl((2.0*n+1)/2/n);
        	f2=sqrtl(2.0*n+1);			
        	 P(i)=f1*cosx*P(j);					// diagonal elements
        	P1(i)=f1*(-sinx*P(j)+P1(j)*cosx);	
         	P(i-1)=f2*sinx*P(j); 					// subdiagonal elements
        	P1(i-1)=f2*(cosx*P(j)+sinx*P1(j));
    	}
    	n=2; m=0; i=2;
    	while(++i<imax-2){
		if(m==n-1){
        		m=0; n++; i++;
        	}
        	else {
        		j=i-n;    					// index for Pn-1,m
        		k=i-2*n+1;					// index for Pn-2,m
            		f3=sqrtl((2.0*n+1)/(n-m)/(n+m));
            		f4=sqrtl(2.0*n-1);
            		f5=sqrtl((n-m-1.0)*(n+m-1.0)/(2.0*n-3));	// remaining elements
             		P(i)=f3*(f4*sinx*P(i-n)-f5*P(k));
            		P1(i)=f3*(f4*(cosx*P(j)+sinx*P1(j))-f5*P1(k));
            		m++;
        	}
    	}
}
void HELP(){
	fprintf(stderr,"\n     The CSHSOFT computes a gravimetric geoid model by using classical Stokes-Helmert method.\n\n");
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"     CSHSOFT -G[model<file>] -A[anomaly<file>] -E[elevation<file>] -T[terrain<file>] ...\n");
	fprintf(stderr,"             ... -M[<value>] -L[<value>] -P[<value>] -R[<value>] -I[<value>] -S -H\n\n");
	fprintf(stderr,"PARAMETERS:\n");
	fprintf(stderr,"     model:     global geopotential model that will be used as a reference model.\n");
	fprintf(stderr,"                It includes harmonic coefficients (n,m,Cnm,Snm,sigmaC,sigmaS).\n\n");
	fprintf(stderr,"     anomaly:   mean free-air gravity anomalies which cover the data area.\n");
	fprintf(stderr,"                It includes grid based data (latitude, longitude, and mean anomaly).\n\n");
	fprintf(stderr,"     elevation: mean topographic elevations which cover the data area.\n");
	fprintf(stderr,"                It includes grid based DTM data (latitude, longitude, and mean elevation).\n\n");
	fprintf(stderr,"     terrain:   terrain corrections which cover the data area.\n");
	fprintf(stderr,"                It includes grid based data (latitude, longitude, and correction).\n\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr,"     -M<value>  maximum expansion of the GGM used in the computation.\n");
	fprintf(stderr,"                default: 630\n\n");
	fprintf(stderr,"     -L<value>  maximum expansion of the Stokes in the computation.\n");
	fprintf(stderr,"                default: 145\n\n");
	fprintf(stderr,"     -P<value>  integration capsize (unit: degree).\n");
	fprintf(stderr,"                default: 0.95\n\n");
	fprintf(stderr,"     -R<value>  limits of the target area (MinLat/MaxLat/MinLon/MaxLon).\n");
	fprintf(stderr,"                default: 45.01/47.00/02.01/4.00\n\n");
	fprintf(stderr,"     -I<value>  intervals of the grid (LatInterval/LonInterval).\n");
	fprintf(stderr,"                default: 0.02/0.02 (72x72 arc-seconds)\n\n");
	fprintf(stderr,"     -S         prints all segments (lat, lon, NGGM, NDg, NH, N), respectively.\n\n");
	fprintf(stderr,"     -H         prints this help.\n\n");
        exit(EXIT_FAILURE);
}
