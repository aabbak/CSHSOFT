/*
 * matris.h, v1.1, is a header file performing various matrix 
 * operations such as inverse, multiplication, addition etc.
 *
 * Copyright (C) 2007 A. Ustun
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
 *           Harita Muh. 42250
 *           Kampus/Konya
 * E-mail  : austun@selcuk.edu.tr
 * 
 * 25.07.2012  Singular vale decomposition functions
 *             svdcmp, svbksb, solvsvd added
 * 22.06.2008  New member functions added
 * 27.02.2007  Creation
 * 
 * Version : 1.1
 * 
 *-----------------------------------------------------------------------------
 */

#ifndef MATRIS_H
#define MATRIS_H


#include <iostream>
#include <stdlib.h>
#include <math.h>

static double maxarg1,maxarg2;
static int iminarg1,iminarg2;
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20;
#define SQR(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
using namespace std;

class matris
{

public:
    double **M;    // 2B isaretci matris
    int n;         // Satir sayisi
    int m;         // Sutun sayisi

    /* Yapilandirici: Vektor *******************************/
    matris(int sat)
    {
        n=sat;
        m=1;
        M=new double*[n];
        if(!M) {cout << "\nVektor yapilandirici: Hafiza ayirma hatasi\n"; exit(EXIT_FAILURE);}
        for(int i=0; i<n; i++)
        {
            M[i]=new double[1];
            M[i][0]=0.0;
        }
    }
    
    /* Yapilandirici: Matris *******************************/
    matris(int sat, int sut)
    {
        n=sat;
        m=sut;
        M = new double*[n];
        if(!M) {cout << "\nMatris yapilandirici: Hafiza ayirma hatasi\n"; exit(EXIT_FAILURE);}
        for(int i=0; i<n; i++)
        {
            M[i]=new double[m];
            for(int j=0;j<m;j++)
                M[i][j]=0.0;
        }
    }
    // Yokedici
    ~matris() {temizlik();}
    
    /* Kopya yapilandirici */
    matris (const matris &mat) {baslat(mat);}
    
    // Matris elemanlarina erisim
    double& operator()(int a, int b) {return *(*(M+a)+b);}
    double& operator()(int a)        {return *(*(M+a)+0);}
    
    /* Atama operatorunun asiri yuklenmesi */
    matris &operator=(const matris &mat)
    {
        if(this!=&mat)
        {
            temizlik();
            baslat(mat);
        }
        return *this;
    }

    /* matris okuma ve yazma */ 
    void oku(const char *);// dosyadan oku
    void yaz(int);         // ekrana  istenen incelikte yaz
    void yaz(char *, int); // dosyaya istenen incelikte yaz


    void random();         // rasgele deger ata
    /*******************************************************/
    /*                  Arkadas Fonksiyonlar               */
    /*******************************************************/
    /* Matris cebri ****************************************/
    friend matris operator+(matris,matris);// Toplama
    friend matris operator-(matris,matris);// Cikarma
    friend matris operator*(matris,matris);// Carpma
    friend matris operator*(double,matris);// Skaler sayi ile carpma
    
    /* Matris islemleri ************************************/
    friend matris transpoze(matris);
    friend double norm2(matris);

    /* Denklem sistemi cozumu *****************************/
    friend void gaussj(matris &,matris &);       // Gauss-Jordan
    friend void lubksb(matris &,int *,matris &); // LU 
    friend void solvlu(matris &,matris &);       // LU ile coz
    friend void solvch(matris,matris &,matris,matris &);// Choloesky
    friend void svdcmp(matris &,matris &,matris &);     // SVD 
    friend void svbksb(matris &,matris &,matris &,matris &,matris &);
    friend matris solvch(matris,matris);
    friend matris vander(matris,matris,matris);
    friend matris solvsvd(matris,matris);

    /* Kare matris islemleri ******************************/
    friend double det(matris);   // LU
    friend double detch(matris); // Cholesky
    friend matris mters(matris); // LU
    friend matris invch(matris); // Cholesky invers

    /* Vektor islemleri ***********************************/
    friend double variance(matris);   // Vektor varyans

    /* Kare matris ayristirma ******************************/
    friend void   ludcmp(matris &, int *, double *); // LU
    friend void   choldc(matris &, double []);       // Cholesky
    friend void   cholsl(matris &, double [], matris &, matris &); // Cholesky yontemine gore denklem cozumu icin ayristirma

    // matris sinifi için cikis operatorunu asiri yukleme
    // << operatoru bir sinifa uye OLAMAZ ancak arkadas olabilir
    friend ostream &operator<<(ostream &, matris);

    // pisagor bagintisi
    friend double pythag(double, double);
private:
/***************************************************************/
    // Nesneyi baslat
    void baslat(const matris &mat)
    {
        n=mat.n;
        m=mat.m;
        M =  new double*[n];
        if(!M) {cout << "baslat: Hafiza ayirma hatasi!!!\n"; exit(EXIT_FAILURE);}
        for(int i=0; i<n; i++)
        {
            M[i]=new double[m];
            pkopya(&M[i][0],&M[i][m-1],&mat.M[i][0]);
        }
    }
    // Nesneyi yoket
    void temizlik()
    {
        for (int i=0;i<n;i++)
        {
            delete [] M[i];
            M[i]=NULL;
        }
        delete [] M;
        M=NULL;
    }
    // Nesneyi belirli bir degerle baslat
    void baslat(double *bptr, double *eptr, double value)
    {
        for(; bptr<=eptr; *bptr++=value)
            ;
    }
    // Isaretci kopyalama
    void pkopya(double *bptr, double *eptr, double *pbcp)
    {
        for(; bptr<=eptr; *bptr++=*pbcp++)
            ;
    }
/***************************************************************/
};

    
#endif

