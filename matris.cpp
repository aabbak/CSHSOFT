/*
 * matris.cpp, v1.0, is a implementation file of matrix class that generates 
 * an object for matrix operations.
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
 * Creation: 20.02.2003
 * Revised : 06.11.2006 A. Ustun    Regularization of the source code 
 *                                  for separate compilation
 * Version : 1.0
 * 
 *-----------------------------------------------------------------------------
 */

#include <iomanip>
#include <fstream>
#include <math.h>
#include "matris.h"
using namespace std;


// Matris elemanlarini dosyadan oku
void matris::oku(const char *dosya)
{
    ifstream dosyaoku(dosya);
    if(dosyaoku.fail()){
        cerr << "oku: " << dosya << " dosyasi bulunamadi!!!\n";exit(EXIT_FAILURE);}
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            dosyaoku >> *(*(M+i)+j);
    dosyaoku.close();
}

// Matris elemanlarini bir dosyaya istenen incelikte yaz
void matris::yaz(char *dosya,int precision)
{
    int i=0;
    int j=0;

    ofstream yazdosya(dosya);

    yazdosya.setf(ios::scientific);
    yazdosya.setf(ios::showpoint);
    yazdosya.precision(precision);

    if(yazdosya.fail()){
        cerr << "yaz: " << dosya << " dosyasi acilamadi!!!";exit(EXIT_FAILURE);}
    if(m<=20)
        for(i=0; i<n; i++)
        {
            for(j=0; j<m; j++)
                yazdosya << setw(8+precision) << *(*(M+i)+j);
            yazdosya << endl;
        }
    else
        for(i=0; i<n; i++)
            for(j=0; j<m; j++)
                yazdosya << setw(8+precision) << *(*(M+i)+j) << endl;
    yazdosya.close();
}

// Matris elemanlarini ekrana yaz
void matris::yaz(int precision)
{
    int i=0;
    int j=0;

    cout.setf(ios::scientific);
    cout.setf(ios::showpoint);
    cout.precision(precision);

    if(m<=20)
        for(i=0; i<n; i++)
        {
            for(j=0; j<m; j++)
                cout << setw(8+precision) << *(*(M+i)+j);
            cout << endl;
        }
    else
        for(i=0; i<n; i++)
            for(j=0; j<m; j++)
                cout << setw(8+precision) << *(*(M+i)+j) << endl;
}

void matris::random()
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
	{
            M[i][j] = rand()/((double)RAND_MAX + 1);
	}
//    double *bptr=&M[0][0];
//    double *eptr=&M[n-1][m-1];
//    for(; bptr<=eptr; *bptr++=Random())
//	    ;
}
// Matris sýnýfý icin cikis operatoru
// myaz ostream tipinde nesne
// mat  matris  tipinde nesne 
// kullaným: cout << A;
// A nesnesi mat nesnesine kopyalaniyor
ostream &operator<<(ostream &myaz, matris mat)
{
    int i=0;
    int j=0;

    // Bilimsel yazým virgulden sonra bes basamak
    myaz.setf(ios::scientific);
    myaz.setf(ios::showpoint);
    myaz.precision(5);
    // sutun sayýsý 20'den az dikdortgen yaz
    if(mat.m<=20)
    {
        for(i=0; i<mat.n; i++)
        {
            for(j=0; j<mat.m; j++)
                // sonuclari 13 karakterlik alana yaz
                myaz << setw(13) << mat(i,j);
            myaz << endl;
        }
    }
    // sutun sayýsý 20'den coksa vektor yaz
    else
    {
        for(i=0; i<mat.n; i++)
            for(j=0; j<mat.m; j++)
                myaz << setw(13) << mat(i,j) << endl;
    }
    return myaz;
}

// Matrisin transpozesini al
matris transpoze(matris A)
{
    int i=0;
    int j=0;

    matris At(A.m, A.n);
    for(i=0; i<At.n; i++)
        for(j=0; j<At.m; j++)
            At(i,j)=A(j,i);
    return At;
}    

// Matrisin normu
double norm2(matris A)
{
    double norm = 0.0;
    for(int i=0;i<A.n;i++)
        for(int j=0;j<A.m;j++)
            norm += A(i,j)*A(i,j);
        
    return sqrt(norm);
}    

// Matris topla
matris operator+(matris A, matris B)
{
    if(A.n!=B.n || A.m!=B.m){
        cerr << "operator+: Matris boyutlari ayni degil!!!\n" << endl;exit(EXIT_FAILURE);}
    for(int i=0;i<A.n;i++)
        for(int j=0;j<A.m;j++)
            A(i,j) += B(i,j);
    return A;
}

// Matris cikar
matris operator-(matris A, matris B)
{
    if(A.n!=B.n || A.m!=B.m){
        cerr << "operator-: Matris boyutlari ayni degil!!!\n" << endl;exit(EXIT_FAILURE);}
    for(int i=0;i<A.n;i++)
        for(int j=0;j<A.m;j++)
            A(i,j) -= B(i,j);
    return A;
}

// Matris carpimi
matris operator*(matris A, matris B)
{
    if(A.m!=B.n){
        cerr << "operator*: Matris boyutlari uygun degil!!!\n" << endl;exit(EXIT_FAILURE);}
    matris C(A.n, B.m);
    for(int i=0;i<A.n;i++)
        for(int j=0;j<B.m;j++)
            for(int k=0;k<A.m;k++)
                C(i,j) += A(i,k)*B(k,j);
    return C;
}

// Skaler sayi ile carpma
matris operator*(double skaler, matris A)
{
    for(int i=0;i<A.n;i++)
        for(int j=0;j<A.m;j++)
            A(i,j) *= skaler;
    return A;
}
double variance(matris v)
{
	double var=0.0;
	for(int i=0;i<v.n;i++)
		var+=v(i)*v(i);
	return var;
}

// Gauss Jordan eleminasyon yontemi ile denklem sistemi cozumu: Ax=B
// Kaynak: Press H, Teukolsky SA, Vetterling WT, Flannery BP (1992) 
// Numerical Recipes in C: The Art of Scientific Computing
// Cambridge University Press, 36--41.
// A(n,n): katsayilar matrisi
// B(n,m): Denklem sisteminin sag taradir. m adet sutundan olusabilir.
//
// YONTEMIN AVANTAJLARI:
// 1. Yontem hem denklem sisteminin cozumunu yapar hem de A matrisinin tersini verir
// 2. Tam pivotlama yaptigi icin sayisal anlamda cozumler cok kararlidir
//
// HANGI DURUMLARDA ONERILMEZ?
// 1. Sadece denklem sistemi cozumu yapilacaksa oteki yontemlere gore uc kat yavas
// 2. Simetrik, pozitif tanimli matrislerin tersi icin CHOLESKY yontemi onerilir
void gaussj(matris &a, matris &b)
{
    int *indxc,*indxr,*ipiv; // pivotlama icin konum belirtecleri
    int i=0,icol=0,irow=0,j=0,k=0,l=0,ll=0;
    double big=.0,dum=.0,pivinv=.0;

    int n=a.n;
    int m=b.m;

    indxc= new int[n];
    indxr= new int[n];
    ipiv = new int[n];

    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) { // indirgenecek sutunlara iliskin ana cevrim
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a(j,k)) >= big) {
                            big=fabs(a(j,k));
                            irow=j;
                            icol=k;
                        }
                    } else if (ipiv[k] > 1) {cerr << "gaussj: Matris tekil-1!!!\n" << endl;exit(EXIT_FAILURE);}
                }
        ++(ipiv[icol]);
        // pivot elemani belirlendi. Gerekli olmasi durumunda kosegene pivot elemanini yerlestirmek icin satirlar arasi
        // degisklik yapilabilir. Sutunlarda degisiklik yapilmaz sadece yeniden isimlendirilir.
        // indxc[i]: i. pivot elemaninin bulundugu sutun
        // indxr[i]: i. pivot elemaninin bulundugu satir
        if (irow != icol) {
            for (l=0;l<n;l++) SWAP(a(irow,l),a(icol,l))
            for (l=0;l<m;l++) SWAP(b(irow,l),b(icol,l))
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a(icol,icol) == 0.0) {cerr << "gaussj: Matris tekil-2!!!\n" << endl;exit(EXIT_FAILURE);}
        pivinv=1.0/a(icol,icol);
        a(icol,icol)=1.0;
        for (l=0;l<n;l++) a(icol,l) *= pivinv;
        for (l=0;l<m;l++) b(icol,l) *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != icol) {
                dum=a(ll,icol);
                a(ll,icol)=0.0;
                for (l=0;l<n;l++) a(ll,l) -= a(icol,l)*dum;
                for (l=0;l<m;l++) b(ll,l) -= b(icol,l)*dum;
            }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                SWAP(a(k,indxr[l]),a(k,indxc[l]));
    }
    delete [] ipiv;
    delete [] indxr;
    delete [] indxc;
}

void ludcmp(matris &a, int *indx, double *d)
{
    int   i=0,imax=0,j=0,k=0;
    double big=.0,dum=.0,sum=.0,temp=.0;
    double *vv;

    int n=a.n;
    vv=new double [n];
    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if ((temp=fabs(a(i,j))) > big) big=temp;
        if (big == 0.0) {cerr << "ludcmp: Matris tekil!!!\n"; exit(EXIT_FAILURE);}
        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            sum=a(i,j);
            for (k=0;k<i;k++) sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            sum=a(i,j);
            for (k=0;k<j;k++)
                sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0;k<n;k++) {
                dum=a(imax,k);
                a(imax,k)=a(j,k);
                a(j,k)=dum;
            }
            *d = -*d;
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a(j,j) == 0.0) a(j,j)=TINY;
        if (j != n-1) {
            dum=1.0/a(j,j);
            for (i=j+1;i<n;i++) a(i,j) *= dum;
        }
    }
    delete [] vv;
}

void lubksb(matris &a, int *indx, matris &b)
//Solves the set of n linear equations A . X = B. Here a[1..n][1..n] is input, not as the matrix
//A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
//as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
//B, and returns with the solution vector X. a, n, and indx are not modied by this routine
//and can be left in place for successive calls with dierent right-hand sides b. This routine takes
//into account the possibility that b will begin with many zero elements, so it is e.cient for use
//in matrix inversion.
{
    int i,ii=0,ip,j;
    double sum;

    int n=a.n;

    for (i=0;i<n;i++) {            //When ii is set to a positive value, it will become the
                                    //index of the rst nonvanishing element of b. We now
                                    //do the forward substitution, equation (2.3.6). The
                                    //only new wrinkle is to unscramble the permutation
                                    //as we go.
        ip=indx[i];
        sum=b(ip);
        b(ip)=b(i);
        if (ii)
            for (j=ii-1;j<=i-1;j++) 
                sum -= a(i,j)*b(j);
        else if (sum)
            ii=i+1;            //A nonzero element was encountered, so from now on we
        b(i)=sum;
    }
    for (i=n-1;i>=0;i--) {            //Now we do the backsubstitution, equation (2.3.7).
        sum=b(i);
        for (j=i+1;j<n;j++) sum -= a(i,j)*b(j);
        b(i)=sum/a(i,i);            //Store a component of the solution vector X.
    }                                //All done!
}

// LUb Denklem sistemini coz
void solvlu(matris &A, matris &b)
{
    if(A.m!=A.n){
        cerr << "solvlu: Kare matris degil!!!\n" << endl;exit(EXIT_FAILURE);}
    int n=A.n;
    int *indx;
    double d;
    indx=new int[n];
    if(!indx){
        cerr << "solvlu: Hafiza ayirma hatasi!!!\n" << endl; exit(EXIT_FAILURE);}
    ludcmp(A,indx,&d);
    lubksb(A,indx,b);
    delete [] indx;
}

matris mters(matris A)
{
    if(A.n!=A.m){
        cerr << "mters: Matris boyutlari uygun degil!!!\n" << endl;exit(EXIT_FAILURE);}
    int i,j;
    int n=A.n;
    int *indx;
    double d;
    indx=new int[n];
    if(!indx){
        cerr << "mters: Hafiza ayirma hatasi!!!\n" << endl; exit(EXIT_FAILURE);}
    matris col(n);
    matris y(n,n);

    ludcmp(A,indx,&d);
    for(j=0; j<n; j++) {
        col(j)=1.0;
        lubksb(A,indx,col);
        for(i=0; i<n; i++) y(i,j)=col(i);
    }
    delete [] indx;
    return y;
}

double det(matris A)
{
    if(A.n!=A.m){
        cerr << "det: Matris boyutlari uygun degil!!!\n" << endl;exit(EXIT_FAILURE);}
    int i;
    int n=A.n;
    int *indx;
    double d;
    indx=new int[n];
    if(!indx){
        cerr << "det: Hafiza ayirma hatasi!!!\n" << endl; exit(EXIT_FAILURE);}
    ludcmp(A,indx,&d);
    for(i=0; i<n; i++)
        d *= A(i,i);
    delete [] indx;
    return d;
}

void choldc(matris &A, double p[])
{
    int i,j,k;
    double sum;
    int n=A.n;
    for(i=0; i<n; i++)
        for(j=i; j<n; j++)
        {
            for (sum=A(i,j), k=i-1; k>=0; k--) sum -= A(i,k)*A(j,k);
            if (i==j) {
                if (sum <= 0.0) {
                    cout << "choldc: A pozitif tanimli matris degil!!!\n"; exit(EXIT_FAILURE);}
                p[i]=sqrt(sum);
            }
            else
                A(j,i) = sum/p[i];
        }
}

void cholsl(matris &A, double p[], matris &b, matris &x)
{
    int i,k;
    double sum;
    int n=A.n;
    for (i=0; i<n; i++)
    {
        for (sum = b(i), k=i-1; k>=0; k--)
            sum -= A(i,k)*x(k);
        x(i) = sum/p[i];
    }
    for (i=n-1; i>=0; i--)
    {
        for (sum = x(i), k=i+1; k<n; k++)
            sum -= A(k,i)*x(k);
        x(i) = sum/p[i];
    }
}

void solvch(matris A, matris &Q, matris b, matris &x)
{
    if (A.n != A.m) {
        cerr << "solvch: Matris kare degil!!!\n"; exit(EXIT_FAILURE);}
    int i,j,k;
    double sum;
    double *p;
    int n=A.n;
    p=new double [n];
    if (!p) {
        cerr << "solvch: Hafiza ayirma hatasi!!!\n"; exit(EXIT_FAILURE);}
    choldc(A,p);
    for(i=n-1; i>=0; i--)
    {
        sum = 1/p[i];
        for(j=i+1; j<n; j++)
            sum -= Q(i,j)*A(j,i);
        Q(i,i) = sum/p[i];
        x(i) += Q(i,i)*b(i);
        for(j=i-1; j>=0; j--)
        {
            sum = 0.0;
            for(k=j+1;k<n;k++)
                sum -= Q(k,i)*A(k,j);
            Q(j,i)=Q(i,j) = sum/p[j];
            x(i) += Q(j,i)*b(j);
            x(j) += Q(i,j)*b(i);
        }
    }
    delete [] p;
}

matris solvch(matris A, matris b)
{
    if (A.n != A.m) {
        cerr << "solvch: Matris kare degil!!!\n"; exit(EXIT_FAILURE);}
    double *p;
    int n=A.n;
    p=new double [n];
    if (!p) {
        cerr << "solvch: Hafiza ayirma hatasi!!!\n"; exit(EXIT_FAILURE);}
    matris x(n);
    choldc(A,p);
    cholsl(A,p,b,x);
    delete [] p;
    return x;
}

double detch(matris A)
{
    if (A.n != A.m) {
        cerr << "detch: Matris kare degil!!!\n"; exit(EXIT_FAILURE);}
    double *p;
    int n=A.n;
    double det=1.0;
    p=new double [n];
    if (!p) {
        cerr << "detch: Hafiza ayirma hatasi!!!\n"; exit(EXIT_FAILURE);}
    choldc(A,p);
    for(int i=0; i<n; i++)
        det *= p[i];
    delete [] p;
    return det*det;
}

matris invch(matris A)
{
    if (A.n != A.m) {
        cerr << "invch: Matris kare degil!!!\n"; exit(EXIT_FAILURE);}
    int i=-1,j=-1,k=-1;
    double sum;
    double *p;
    int n=A.n;
    p=new double [n];
    if (!p) {
        cerr << "invch: Hafiza ayirma hatasi!!!\n"; exit(EXIT_FAILURE);}
    choldc(A,p);
    matris Q(n,n);
    for(i=n-1; i>=0; i--)
    {
        sum = 1/p[i];
        for(j=i+1; j<n; j++)
            sum -= Q(i,j)*A(j,i);
        Q(i,i) = sum/p[i];
        for(j=i-1; j>=0; j--)
        {
            sum = 0.0;
            for(k=j+1;k<n;k++)
                sum -= Q(k,i)*A(k,j);
            Q(j,i)=Q(i,j) = sum/p[j];
        }
    }
    delete [] p;
    return Q;
}

matris vander(matris x, matris q)
{
    int i,j,k;
    int n=x.n;
    double xx,b,s,t;

    matris w(n);
    matris c(n);

    if (n==1) w(0)=q(0);
    else {
        c(n-1) = -x(0);
        for(i=1; i<n; i++) {
            xx = -x(i);
            for (j=(n-1-i); j<(n-1); j++) c(j) += xx*c(j+1);
            c(n-1) += xx;
        }
        for(i=0; i<n; i++) {
            xx=x(i);
            t=b=1.0;
            s=q(n-1);
            for (k=n-1; k>=1; k--) {
                b  = c(k)+xx*b;
                s += q(k-1)*b;
                t  = xx*t+b;
            }
            w(i) = s/t;
        }
    }
    return w;
}

void svdcmp(matris &a, matris &w, matris &v)
{
    int n=a.m;  // sutun
    int m=a.n;  // satir
	double pythag(double a, double b);
	int flag=0,i=0,its=0,j=0,jj=0,k=0,l=0,nm=0;
	double anorm=0.,c=0.,f=0.,g=0.,h=0.,s=0.,scale=0.,x=0.,y=0.,z=0.;
    
    matris rv1(n);
	//rv1=vector(1,n);

	g=scale=anorm=0.0;
	//for (i=1;i<=n;i++) {
	for (i=0;i<n;i++) {
		l=i+1;
		rv1(i)=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a(k,i));
			if (scale) {
				for (k=i;k<m;k++) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f=a(i,i);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,i)=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
					f=s/h;
					for (k=i;k<m;k++) a(k,j) += f*a(k,i);
				}
				for (k=i;k<m;k++) a(k,i) *= scale;
			}
		}
		w(i)=scale *g;
		g=s=scale=0.0;
		if (i < m && i != (n-1)) {
			for (k=l;k<n;k++) scale += fabs(a(i,k));
			if (scale) {
				for (k=l;k<n;k++) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f=a(i,l);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,l)=f-g;
				for (k=l;k<n;k++) rv1(k)=a(i,k)/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(j,k)*a(i,k);
					for (k=l;k<n;k++) a(j,k) += s*rv1(k);
				}
				for (k=l;k<n;k++) a(i,k) *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w(i))+fabs(rv1(i))));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<n;j++)
					v(j,i)=(a(i,j)/a(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1(i);
		l=i;
	}
	for (i=(IMIN(m,n)-1);i>=0;i--) {
		l=i+1;
		g=w(i);
		for (j=l;j<n;j++) a(i,j)=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
				f=(s/a(i,i))*g;
				for (k=i;k<m;k++) a(k,j) += f*a(k,i);
			}
			for (j=i;j<m;j++) a(j,i) *= g;
		} else for (j=i;j<m;j++) a(j,i)=0.0;
		a(i,i)=1+a(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(fabs(rv1(l))+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w(nm))+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1(i);
					rv1(i)=c*rv1(i);
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w(i);
					h=pythag(f,g);
					w(i)=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a(j,nm);
						z=a(j,i);
						a(j,nm)=y*c+z*s;
						a(j,i)=z*c-y*s;
					}
				}
			}
			z=w(k);
			if (l == k) {
				if (z < 0.0) {
					w(k) = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == 30) fprintf(stderr,"no convergence in 30 svdcmp iterations\n");
			x=w(l);
			nm=k-1;
			y=w(nm);
			g=rv1(nm);
			h=rv1(k);
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1(i);
				y=w(i);
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1(j)=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=pythag(f,h);
				w(j)=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a(jj,j);
					z=a(jj,i);
					a(jj,j)=y*c+z*s;
					a(jj,i)=z*c-y*s;
				}
			}
			rv1(l)=0.0;
			rv1(k)=f;
			w(k)=x;
		}
	}
	//free_vector(rv1,1,n);
}

void svbksb(matris &u, matris &w, matris &v, matris &b, matris &x)
{
    int n=u.m;  // sutun
    int m=u.n;  // satir
    int jj=0,j=0,i=0;
    double s=0.;
    
    matris tmp(n);
    
    for (j=0;j<n;j++) {
        s=0.0;
        if (w(j)) {
            for (i=0;i<m;i++) s += u(i,j)*b(i);
            s /= w(j);
        }
        tmp(j)=s;
    }
    for (j=0;j<n;j++) {
        s=0.0;
        for (jj=0;jj<n;jj++) s += v(j,jj)*tmp(jj);
        x(j)=s;
    }
}
matris solvsvd(matris U,matris b)
{
    int n=U.m;
    int j=0;
    matris w(n);
    matris V(n,n);
    matris x(n);

    svdcmp(U,w,V);

    double wmax=0.0;
    for(j=0;j<n;j++) if (w(j)>wmax) wmax=w(j);
    double wmin=wmax*1.0e-15;
    for(j=0;j<n;j++) if (w(j)<wmin) w(j)=0.0;

    svbksb(U,w,V,b,x);

    return x;
}

double pythag(double a, double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
