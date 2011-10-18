/* [response orientation nms filterBank] = steerableDetector(image, filterOrder, sigma);
 *
 * (c) Francois Aguet, 07/12/2011 (last modified Jul 14, 2011). Adapted from SteerableJ package, 2008.
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a steerableDetector.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.14" -I"..\..\mex\include" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output steerableDetector steerableDetector.cpp
 */

#include <math.h>
#include <string.h>
#include <gsl/gsl_poly.h>

#include "mex.h"
#include "conv2D.h"

#define PI 3.141592653589793
#define TOL 1e-12


static double approxZero(double n) {
    if (fabs(n) < TOL) {
        return 0.0;
    } else {
        return n;
    }
}

static double opposite(double theta) {
    if (theta > 0.0) { // values in (-PI,PI]
        return theta - PI;
    } else {
        return theta + PI;
    }
}

static int getTemplateN(int M) {
    if (M==1) {
        return M*(M+3)/2;
    } else {
        return M*(M+3)/2 - getTemplateN(M-1);
    }
}


double* getWeights(int M, double sigma) {
    double* alpha;
    switch (M) {
        case 1:
            alpha = new double[1];
            alpha[0] = -sqrt(2.0/PI);
            break;
        case 2:
            alpha = new double[2];
            alpha[0] = 0.16286750396763996;
            alpha[1] = -0.4886025119029199;
            break;
        case 3: // mu = 0
            alpha = new double[3];
            alpha[0] = -0.966;
            alpha[1] = -0.256;
            alpha[2] = 0.0; 
            break;
        case 4: // mu = 0.25
            alpha = new double[5];
            alpha[0] = 0.113;
            alpha[1] = -0.392;
            alpha[2] = 0.025;
            alpha[3] = -0.184;
            alpha[4] = 0.034; 
            break;
        case 5:
            alpha = new double[5];
            alpha[0] = -1.1215;
            alpha[1] = -0.5576;
            alpha[2] = -0.018;
            alpha[3] = -0.0415;
            alpha[4] = -0.0038;
            break;
        default:
            alpha = NULL;
    }
    return alpha;
}


void computeBaseTemplates(double* input, int nx, int ny, int M, double sigma, double** templates) {
    
    int wWidth = (int)(4.0*sigma);
    int nk = wWidth+1;
    
    double* aKernel = new double[nk];
    double* bKernel = new double[nk];
    double* g = new double[nk];
    double* buffer = new double[nx*ny];
    double d;
    double sigma2 = sigma*sigma;
    double sigma4 = sigma2*sigma2;
    double sigma6 = sigma4*sigma2;
    double sigma8 = sigma4*sigma4;
    double sigma10 = sigma6*sigma4;
    double sigma12 = sigma8*sigma4;
    
    for (int i=0;i<nk;i++) {
        g[i] = exp(-(i*i)/(2.0*sigma2));
    }
    
    if (M == 1 || M == 3 || M == 5) {
        
        d = 2.0*PI*sigma4;
        for (int i=0;i<nk;i++) {
            aKernel[i] = -i*g[i] / d;
        }
        
        // g_x
        convolveOddX(input, nx, ny, aKernel, nk, buffer);
        convolveEvenY(buffer, nx, ny, g, nk, templates[0]);
        
        // g_y
        convolveOddY(input, nx, ny, aKernel, nk, buffer);
        convolveEvenX(buffer, nx, ny, g, nk, templates[1]);
        
        if (M == 3 || M == 5) {
            
            d = 2.0*PI*sigma8;
            for (int i=0;i<nk;i++) {
                aKernel[i] = (3.0*i*sigma2 - i*i*i) * g[i] / d;
            }
            // g_xxx
            convolveOddX(input, nx, ny, aKernel, nk, buffer);
            convolveEvenY(buffer, nx, ny, g, nk, templates[2]);
            // g_yyy
            convolveOddY(input, nx, ny, aKernel, nk, buffer);
            convolveEvenX(buffer, nx, ny, g, nk, templates[5]);
            for (int i=0;i<nk;i++) {
                aKernel[i] = (sigma2 - i*i) * g[i] / d;
                bKernel[i] = i*g[i];
            }
            // gxxy
            convolveEvenX(input, nx, ny, aKernel, nk, buffer);
            convolveOddY(buffer, nx, ny, bKernel, nk, templates[3]);
            // gxyy
            convolveEvenY(input, nx, ny, aKernel, nk, buffer);
            convolveOddX(buffer, nx, ny, bKernel, nk, templates[4]);
        }
        if (M == 5) {
            
            d = 2.0*PI*sigma12;
            for (int i=0;i<nk;i++) {
                aKernel[i] = -i*(i*i*i*i - 10.0*i*i*sigma2 + 15.0*sigma4) * g[i] / d;
            }
            // gxxxxx
            convolveOddX(input, nx, ny, aKernel, nk, buffer);
            convolveEvenY(buffer, nx, ny, g, nk, templates[6]);
            // gyyyyy
            convolveOddY(input, nx, ny, aKernel, nk, buffer);
            convolveEvenX(buffer, nx, ny, g, nk, templates[11]);
            for (int i=0;i<nk;i++) {
                aKernel[i] = (i*i*i*i - 6.0*i*i*sigma2 + 3.0*sigma4) * g[i] / d;
                bKernel[i] = -i * g[i];
            }
            // g_xxxxy
            convolveEvenX(input, nx, ny, aKernel, nk, buffer);
            convolveOddY(buffer, nx, ny, bKernel, nk, templates[7]);
            // g_xyyyy
            convolveEvenY(input, nx, ny, aKernel, nk, buffer);
            convolveOddX(buffer, nx, ny, bKernel, nk, templates[10]);
            for (int i=0;i<nk;i++) {
                aKernel[i] = i*(i*i - 3.0*sigma2) * g[i] / d;
                bKernel[i] = (sigma2 - i*i) * g[i];
            }
            // g_xxxyy
            convolveOddX(input, nx, ny, aKernel, nk, buffer);
            convolveEvenY(buffer, nx, ny, bKernel, nk, templates[8]);
            // g_xxyyy
            convolveOddY(input, nx, ny, aKernel, nk, buffer);
            convolveEvenX(buffer, nx, ny, bKernel, nk, templates[9]);
        }
    } else { //(M == 2 || M == 4)
        
        d = 2.0*PI*sigma6;
        for (int i=0;i<nk;i++) {
            aKernel[i] = (i*i - sigma2) * g[i] / d;
        }
        // g_xx
        convolveEvenX(input, nx, ny, aKernel, nk, buffer);
        convolveEvenY(buffer, nx, ny, g, nk, templates[0]);
        // g_yy
        convolveEvenY(input, nx, ny, aKernel, nk, buffer);
        convolveEvenX(buffer, nx, ny, g, nk, templates[2]);
        for (int i=0;i<nk;i++) {
            aKernel[i] = i * g[i];
            bKernel[i] = aKernel[i] / d;
        }
        // g_xy
        convolveOddX(input, nx, ny, aKernel, nk, buffer);
        convolveOddY(buffer, nx, ny, bKernel, nk, templates[1]);
        
        if (M == 4) {
            
            d = 2.0*PI*sigma10;
            for (int i=0;i<nk;i++) {
                aKernel[i] = (i*i*i*i - 6.0*i*i*sigma2 + 3.0*sigma4) * g[i] / d;
            }
            // g_xxxx
            convolveEvenX(input, nx, ny, aKernel, nk, buffer);
            convolveEvenY(buffer, nx, ny, g, nk, templates[3]);
            // g_yyyy
            convolveEvenY(input, nx, ny, aKernel, nk, buffer);
            convolveEvenX(buffer, nx, ny, g, nk, templates[7]);
            for (int i=0;i<nk;i++) {
                aKernel[i] = i * (i*i - 3.0*sigma2) * g[i] / d;
                bKernel[i] = i * g[i];
            }
            // g_xxxy
            convolveOddX(input, nx, ny, aKernel, nk, buffer);
            convolveOddY(buffer, nx, ny, bKernel, nk, templates[4]);
            // g_xyyy
            convolveOddY(input, nx, ny, aKernel, nk, buffer);
            convolveOddX(buffer, nx, ny, bKernel, nk, templates[6]);
            for (int i=0;i<nk;i++) {
                aKernel[i] = (sigma2 - i*i) * g[i];
                bKernel[i] = aKernel[i] / d;
            }
            // g_xxyy
            convolveEvenX(input, nx, ny, aKernel, nk, buffer);
            convolveEvenY(buffer, nx, ny, bKernel, nk, templates[5]);
        }
    }
    
    // free memory
    delete[] aKernel;
    delete[] bKernel;
    delete[] g;
    delete[] buffer;
}




// Point responses
/*double pointRespM1(int i, double angle, double* alpha, double** templates) {
    
    double cosT = cos(angle);
    double sinT = sin(angle);
    
    double gxi = templates[0][i];
    double gyi = templates[1][i];
    double a11 = alpha[0];
    
    double result = a11 * (cosT*gxi + sinT*gyi);
    return result;
}*/


double pointRespM2(int i, double angle, double* alpha, double** templates) {
    
    double cosT = cos(angle);
    double sinT = sin(angle);
    
    double gxxi = templates[0][i];
    double gxyi = templates[1][i];
    double gyyi = templates[2][i];
    double a20 = alpha[0];
    double a22 = alpha[1];
    
    double result = cosT*cosT*(a22*gyyi + a20*gxxi)
                  + cosT*sinT*2.0*(a20-a22)*gxyi
                  + sinT*sinT*(a22*gxxi + a20*gyyi);
    return result;
}


double pointRespM3(int i, double angle, double* alpha, double** templates) {
    
    double cosT = cos(angle);
    double sinT = sin(angle);
    double cosT2 = cosT*cosT;
    double sinT2 = sinT*sinT;
    
    double gxi = templates[0][i];
    double gyi = templates[1][i];
    double gxxxi = templates[2][i];
    double gxxyi = templates[3][i];
    double gxyyi = templates[4][i];
    double gyyyi = templates[5][i];
    
    double a10 = alpha[0];
    double a30 = alpha[1];
    double a32 = alpha[2];
    
    double result = a10*(sinT*gyi + cosT*gxi)
            + sinT2*sinT*(a30*gyyyi + a32*gxxyi)
            + cosT*sinT2*(a32*gxxxi + (3.0*a30-2.0*a32)*gxyyi)
            + cosT2*sinT*(a32*gyyyi + (3.0*a30-2.0*a32)*gxxyi)
            + cosT2*cosT*(a32*gxyyi + a30*gxxxi);
    
    return result;
}


double pointRespM4(int i, double angle, double* alpha, double** templates) {
    
    double cosT = cos(angle);
    double sinT = sin(angle);
    double cosTsinT = cosT*sinT;
    double cosT2 = cosT*cosT;
    double sinT2 = sinT*sinT;
    
    double gxxi = templates[0][i];
    double gxyi = templates[1][i];
    double gyyi = templates[2][i];    
    double gxxxxi = templates[3][i];
    double gxxxyi = templates[4][i];
    double gxxyyi = templates[5][i];
    double gxyyyi = templates[6][i];
    double gyyyyi = templates[7][i];
    
    double a20 = alpha[0];
    double a22 = alpha[1];
    double a40 = alpha[2];
    double a42 = alpha[3];
    double a44 = alpha[4];            
    
    double a = (a20-a22)*gxyi;
    
    double result = cosT2*cosT2 * (a20*gxxi+a22*gyyi+a40*gxxxxi+a42*gxxyyi+a44*gyyyyi)
    + cosT2*cosTsinT * 2.0*(a+2.0*a40*gxxxyi+a42*(gxyyyi-gxxxyi)-2.0*a44*gxyyyi)
    + cosT2*sinT2 * ((a20+a22)*(gxxi+gyyi)+a42*(gxxxxi+gyyyyi)+(6.0*(a40+a44)-4.0*a42)*gxxyyi)
    + sinT2*cosTsinT * 2.0*(a+2.0*a40*gxyyyi+a42*(gxxxyi-gxyyyi)-2.0*a44*gxxxyi)
    + sinT2*sinT2 * (a20*gyyi+a22*gxxi+a40*gyyyyi+a42*gxxyyi+a44*gxxxxi);
    return result;
}


double pointRespM5(int i, double angle, double* alpha, double** templates) {
    
    double cosT = cos(angle);
    double sinT = sin(angle);
    double cosT2 = cosT*cosT;
    double sinT2 = sinT*sinT;
    double cosT3 = cosT2*cosT;
    double sinT3 = sinT2*sinT;
    double sinT4 = sinT2*sinT2;
    double cosT4 = cosT2*cosT2;
    double sinT5 = sinT2*sinT3;
    double cosT5 = cosT2*cosT3;
    
    double gxi = templates[0][i];
    double gyi = templates[1][i];
    double gxxxi = templates[2][i];
    double gxxyi = templates[3][i];
    double gxyyi = templates[4][i];
    double gyyyi = templates[5][i];
    double gxxxxxi = templates[6][i];
    double gxxxxyi = templates[7][i];
    double gxxxyyi = templates[8][i];
    double gxxyyyi = templates[9][i];
    double gxyyyyi = templates[10][i];
    double gyyyyyi = templates[11][i];
    
    double a10 = alpha[0];
    double a30 = alpha[1];
    double a32 = alpha[2];
    double a52 = alpha[3];
    double a54 = alpha[4];
    
    /*double result = cosT2*cosT3 * (a11*gyi + a31*gxxyi + a33*gyyyi + a51*gxxxxyi + a53*gxxyyyi)
    + cosT2*cosT2*sinT * (-a11*gxi - a31*gxxxi - (3.0*a33-2.0*a31)*gxyyi + a51*(4.0*gxxxyyi - gxxxxxi) + a53*(2.0*gxyyyyi - 3.0*gxxxyyi))
    + cosT3*sinT2 * ( 2.0*a11*gyi + (3.0*a33-a31)*gxxyi + (a33+a31)*gyyyi + a51*(6.0*gxxyyyi-4.0*gxxxxyi) + a53*(gyyyyyi-6.0*gxxyyyi+3.0*gxxxxyi))
    + cosT2*sinT3 * (-2.0*a11*gxi - (3.0*a33-a31)*gxyyi - (a33+a31)*gxxxi - a51*(6.0*gxxxyyi-4.0*gxyyyyi) - a53*(gxxxxxi-6.0*gxxxyyi+3.0*gxyyyyi))
    + cosT*sinT2*sinT2 * ( a11*gyi + a31*gyyyi + (3.0*a33-2.0*a31)*gxxyi - a51*(4.0*gxxyyyi - gyyyyyi) - a53*(2.0*gxxxxyi - 3.0*gxxyyyi))
    - sinT2*sinT3 * (a11*gxi + a31*gxyyi + a33*gxxxi + a51*gxyyyyi + a53*gxxxyyi);*/
    
    /*double result = a10*sinT*gyi + a30*sinT3*gyyyi + a32*sinT3*gxxyi
            + a52*sinT5*gxxyyyi + cosT5*(a54*gxyyyyi + a52*gxxxyyi) + a54*sinT5*gxxxxyi
            + cosT4*sinT*(a54*gyyyyyi + 3.0*a52*gxxyyyi - 4.0*a54*gxxyyyi - 2.0*a52*gxxxxyi)
            + cosT2*(a32*sinT*gyyyi + a52*sinT3*gyyyyyi + 3.0*a30*sinT*gxxyi - 2.0*a32*sinT*gxxyi
            - 6.0*a52*sinT3*gxxyyyi + 6.0*a54*sinT3*gxxyyyi + 3.0*a52*sinT3*gxxxxyi - 4.0*a54*sinT3*gxxxxyi)
            + cosT3*(a32*gxyyi + 3.0*a52*sinT2*gxyyyyi - 4.0*a54*sinT2*gxyyyyi + a30*gxxxi
            - 6.0*a52*sinT2*gxxxyyi + 6.0*a54*sinT2*gxxxyyi + a52*sinT2*gxxxxxi)
            + cosT*(a10*gxi + 3.0*a30*sinT2*gxyyi - 2.0*a32*sinT2*gxyyi - 2.0*a52*sinT4*gxyyyyi
            + a32*sinT2*gxxxi + 3.0*a52*sinT4*gxxxyyi - 4.0*a54*sinT4*gxxxyyi + a54*sinT4*gxxxxxi);*/
    
    double result = a10*sinT*gyi
            + sinT3*(a30*gyyyi + a32*gxxyi)
            + sinT5*(a52*gxxyyyi + a54*gxxxxyi)
            + cosT5*(a54*gxyyyyi + a52*gxxxyyi) 
            + cosT4*sinT*(a54*gyyyyyi + 3.0*a52*gxxyyyi - 4.0*a54*gxxyyyi - 2.0*a52*gxxxxyi)
            
            + cosT2*(a32*sinT*gyyyi + a52*sinT3*gyyyyyi + 3.0*a30*sinT*gxxyi - 2.0*a32*sinT*gxxyi
            - 6.0*a52*sinT3*gxxyyyi + 6.0*a54*sinT3*gxxyyyi + 3.0*a52*sinT3*gxxxxyi - 4.0*a54*sinT3*gxxxxyi)
            + cosT3*(a32*gxyyi + 3.0*a52*sinT2*gxyyyyi - 4.0*a54*sinT2*gxyyyyi + a30*gxxxi
            - 6.0*a52*sinT2*gxxxyyi + 6.0*a54*sinT2*gxxxyyi + a52*sinT2*gxxxxxi)
            + cosT*(a10*gxi + 3.0*a30*sinT2*gxyyi - 2.0*a32*sinT2*gxyyi - 2.0*a52*sinT4*gxyyyyi
            + a32*sinT2*gxxxi + 3.0*a52*sinT4*gxxxyyi - 4.0*a54*sinT4*gxxxyyi + a54*sinT4*gxxxxxi);

    return result;
}




int getRealRoots(double* z, int nz, double* roots) {
    int nr = nz/2; // total roots

    int nrr = 0; // real roots
    for (int k=0;k<nr;++k) {
        if (z[2*k+1]==0.0)
            nrr++;
    }
    roots = new double[nrr];
    nrr = 0;
    for (int k=0;k<nr;++k) {
        if (z[2*k+1]==0.0)
            roots[nrr++] = z[2*k];
    }
    return nrr;
}



void filterM1(double** templates, int nx, int ny, double* alpha, double* response, double* orientation) {
    
    double* gx = templates[0];
    double* gy = templates[1];
    double a11 = alpha[0];
    
    double* tRoots = new double[2];
    double gxi, gyi;
    double temp;
    
    for (int i=0;i<nx*ny;++i) {
        gxi = approxZero(gx[i]);
        gyi = approxZero(gy[i]);
        
        if ( gxi==0.0 && gyi==0.0 ) {
            response[i] = 0.0;
            orientation[i] = 0.0;
        } else {
            if (gyi == 0.0) {
                tRoots[0] = 0.0;
                tRoots[1] = PI;
            } else {
                tRoots[0] = atan(gyi/gxi);
                tRoots[1] = opposite(tRoots[0]);
            }
            orientation[i] = tRoots[0];
            response[i] = a11 * (cos(tRoots[0])*gxi + sin(tRoots[0])*gyi);
            temp = a11 * (cos(tRoots[1])*gxi + sin(tRoots[1])*gyi);
            if (temp > response[i]) {
                response[i] = temp;
                orientation[i] = tRoots[1];
            }
        }
    }
    delete[] tRoots;
}




// quadratic root solution
void filterM2(double** templates, int nx, int ny, double* alpha, double* response, double* orientation) {
    
    double* gxx = templates[0];
    double* gxy = templates[1];
    double* gyy = templates[2];
    double a20 = alpha[0];
    double a22 = alpha[1];
    
    double A, B, C;
    double a = a20-a22;
    double temp;

    for (int i=0;i<nx*ny;++i) {
        
        C = a*gxy[i];
        B = a*(gyy[i]-gxx[i]);
        A = -C;
        
        if (A == 0.0) { // -> linear
            if (B == 0.0) { // -> null, solve
                orientation[i] = 0.0;
                response[i] = pointRespM2(i, 0.0, alpha, templates);
            } else { // solve linear
                if (C == 0.0) {
                    orientation[i] = 0.0;
                    response[i] = pointRespM2(i, 0.0, alpha, templates);
                    temp = pointRespM2(i, PI/2.0, alpha, templates);
                    if (temp > response[i]) {
                        response[i] = temp;
                        orientation[i] = PI/2.0;
                    }
                } else {
                    orientation[i] = atan(-C/B);
                    response[i] = pointRespM2(i, orientation[i], alpha, templates);
                    temp = pointRespM2(i, orientation[i]+PI/2.0, alpha, templates);
                    if (temp > response[i]) {
                        response[i] = temp;
                        orientation[i] += PI/2.0;
                    }
                }
            }
        } else { // solve quadratic
            double* xRoots = new double[2];
            gsl_poly_solve_quadratic (A, B, C, &xRoots[0], &xRoots[1]);
            
            double* tRoots = new double[2];
            tRoots[0] = atan(xRoots[0]);
            tRoots[1] = atan(xRoots[1]);
            response[i] = pointRespM2(i, tRoots[0], alpha, templates);
            orientation[i] = tRoots[0];
            temp = pointRespM2(i, tRoots[1], alpha, templates);
            if (temp > response[i]) {
                response[i] = temp;
                orientation[i] = tRoots[1];
            }
            delete[] xRoots;
            delete[] tRoots;
        }
    }
}



void filterM3(double** templates, int nx, int ny, double* alpha, double* response, double* orientation) {
    
    double* gx = templates[0];
    double* gy = templates[1];
    double* gxxx = templates[2];
    double* gxxy = templates[3];
    double* gxyy = templates[4];
    double* gyyy = templates[5];
    
    /*double a11 = alpha[0];
    double a31 = alpha[1];
    double a33 = alpha[2];*/
    
    double a10 = alpha[0];
    double a30 = alpha[1];
    double a32 = alpha[2];
    
    double A, B, C, D;
    
    //double a1 = 2.0*a31-3.0*a33;
    //double a2 = 6.0*a33-7.0*a31;
    //double g1, g2;
    
    int nr, nt, deg;
    
    // Derivative of the steering function           
    /*sinT3 * (-a10*gxi - 3.0*a30*gxyyi + 2.0*a32*gxyyi - a32*gxxxi)
    + sinT2*cosT * (a10*gyi + 3.0*a30*gyyyi - 2.0*a32*gyyyi - 6.0*a30*gxxyi + 7.0*a32*gxxyi)
    + cosT2*sinT * (-a10*gxi + 6.0*a32*gxyyi - 7.0*a32*gxyyi - 2.0*a30*gxxxi + 2.0*a32*gxxxi)
    + cosT3 * (a10*gyi + a32*gyyyi*3.0*a30*gxxyi-3.0*a32*gxxyi);*/
    
    
    
    for (int i=0;i<nx*ny;++i) {
        
        /*g1 = -a11*gy[i];
        g2 = -a11*gx[i];
        
        A = g1 - a31*gyyy[i] + a1*gxxy[i];
        B = g2 + a2*gxyy[i] + a1*gxxx[i];
        C = g1 + a2*gxxy[i] + a1*gyyy[i];
        D = g2 - a31*gxxx[i] + a1*gxyy[i];*/
        
        A = -a10*gx[i] + (2.0*a32- 3.0*a30)*gxyy[i] - a32*gxxx[i];
        B =  a10*gy[i] + (3.0*a30-2.0*a32)*gyyy[i] - (6.0*a30-7.0*a32)*gxxy[i];
        C = -a10*gx[i] + (2.0*a32-3.0*a30)*gxxx[i] + (6.0*a30-7.0*a32)*gxyy[i]; 
        D =  a10*gy[i] + a32*gyyy[i] + (3.0*a30-2.0*a32)*gxxy[i];
        
        A = approxZero(A);
        B = approxZero(B);
        C = approxZero(C);
        D = approxZero(D);
        
        double* roots;
        double* z;
        if (A == 0.0) { // -> quadratic
            if (B == 0.0) { // -> linear
                if (C == 0.0) {// -> null, fixed solution
                    deg = 1;
                    z = new double[2*deg];
                    z[0] = 0.0;
                    z[1] = 0.0;
                } else {
                    deg = 1;
                    double a[2] = {D, C};
                    z = new double[2*deg];
                    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                    gsl_poly_complex_solve(a, deg+1, w, z);
                    gsl_poly_complex_workspace_free(w);
                }
            } else { // B!=0
                deg = 2;
                double a[3] = {D, C, B};
                z = new double[2*deg];
                gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                gsl_poly_complex_solve(a, deg+1, w, z);
                gsl_poly_complex_workspace_free(w);
            }
        } else { // solve cubic
            deg = 3;
            double a[4] = {D, C, B, A};
            z = new double[2*deg];
            gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
            gsl_poly_complex_solve(a, deg+1, w, z);
            gsl_poly_complex_workspace_free(w);
        }
    
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0)
                nr++;
        }
        roots = new double[nr];
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0) {
                roots[nr++] = z[2*k];
            }
        }
        //nr = getRealRoots(z, 2*deg, roots);
        
        double* tRoots;
        if (nr == 0) {
            nt = 4;
            tRoots = new double[nt];
            tRoots[0] = -PI/2.0;
            tRoots[1] = 0.0;
            tRoots[2] = PI/2.0;
            tRoots[3] = PI;
        } else {
            nt = 2*nr;
            tRoots = new double[nt];
            for (int k=0;k<nr;k++) {
                tRoots[k] = atan(roots[k]);
                tRoots[k+nr] = opposite(tRoots[k]);
            }
        }
        delete[] roots;
        
        
        response[i] = pointRespM3(i, tRoots[0], alpha, templates);
        orientation[i] = tRoots[0];
        
        double temp;
        for (int k=1;k<nt;k++) {
            temp = pointRespM3(i, tRoots[k], alpha, templates);
            if (temp > response[i]) {
                response[i] = temp;
                orientation[i] = tRoots[k];
            }
        }
        delete[] tRoots;
        delete[] z;
    }
}




void filterM4(double** templates, int nx, int ny, double* alpha, double* response, double* orientation) {
    
    double* gxx = templates[0];
    double* gxy = templates[1];
    double* gyy = templates[2];
    double* gxxxx = templates[3];
    double* gxxxy = templates[4];
    double* gxxyy = templates[5];
    double* gxyyy = templates[6];
    double* gyyyy = templates[7];

    double a20 = alpha[0];
    double a22 = alpha[1];
    double a40 = alpha[2];
    double a42 = alpha[3];
    double a44 = alpha[4];
    
    double g1, g2, g3;
    double A, B, C, D, E;
    
    double a1 = 2.0*a44-a42;
    double a2 = 2.0*a40-a42;
    double a3 = a22-a20;
    double a4 = 6.0*(a44-a42+a40);
    
    int nr, deg;
    double delta;
    
    for (int i=0;i<nx*ny;i++) {
        
        g1 = a3*gxy[i];
        g2 = a3*(gxx[i]-gyy[i]);
        g3 = a4*gxxyy[i];
        
        A = a1*gxxxy[i] - a2*gxyyy[i] + g1;
        B = a1*gxxxx[i] + a2*gyyyy[i] - g3 + g2;
        C = a4*(gxyyy[i] - gxxxy[i]);
        D = -a1*gyyyy[i] - a2*gxxxx[i] + g3 + g2;
        E = a2*gxxxy[i] - a1*gxyyy[i] - g1;
        
        A = approxZero(A);
        C = approxZero(C);
        E = approxZero(E);
        
        double* roots;
        double* z;
        if (A == 0.0) { // -> cubic
            if (B == 0.0) { // -> quadratic
                if (C == 0.0) { // -> linear
                    if (D == 0.0) { // solve null
                        deg = 1;
                        z = new double[2*deg];
                        z[0] = 0.0;
                        z[1] = 0.0;
                    } else { // solve linear
                        deg = 1;
                        double a[2] = {E, D};
                        z = new double[2*deg];
                        gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                        gsl_poly_complex_solve(a, deg+1, w, z);
                        gsl_poly_complex_workspace_free(w);
                    }
                } else { // solve quadratic
                    deg = 2;
                    double a[3] = {E, D, C};
                    z = new double[2*deg];
                    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                    gsl_poly_complex_solve(a, deg+1, w, z);
                    gsl_poly_complex_workspace_free(w);
                }
            } else { // solve cubic
                if ( (C == 0.0) && (E == 0.0) ) {
                    delta = -D/B;
                    if (delta > 0.0) {
                        deg = 3;
                        delta = sqrt(delta);
                        z = new double[2*deg];
                        z[0] = 0.0;
                        z[1] = 0.0;
                        z[2] = delta;
                        z[3] = 0.0;
                        z[4] = -delta;
                        z[5] = 0.0;
                    } else {
                        deg = 1;
                        z = new double[2*deg];
                        z[0] = 0.0;
                        z[1] = 0.0;
                    }
                } else {
                    deg = 3;
                    double a[4] = {E, D, C, B};
                    z = new double[2*deg];
                    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                    gsl_poly_complex_solve(a, deg+1, w, z);
                    gsl_poly_complex_workspace_free(w);
                }
            }
        } else { // solve quartic
            deg = 4;
            double a[5] = {E, D, C, B, A};
            z = new double[2*deg];
            gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
            gsl_poly_complex_solve(a, deg+1, w, z);
            gsl_poly_complex_workspace_free(w);
        }
            
        // # real roots
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0)
                nr++;
        }
        roots = new double[nr];
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0) {
                roots[nr++] = z[2*k];
            }
        }

        
        double* tRoots;
        if (nr == 0) { //orientation[i] = 0.0;
            nr = 2;
            tRoots = new double[nr];
            tRoots[0] = 0.0;
            tRoots[1] = PI/2.0;
        } else {
            if (roots[0] == 0.0) {
                nr++;
                tRoots = new double[nr];
                for (int k=0;k<nr-1;k++) {
                    tRoots[k] = atan(roots[k]);
                }
                tRoots[nr-1] = PI/2;
            } else {
                tRoots = new double[nr];
                for (int k=0;k<nr;k++) {
                    tRoots[k] = atan(roots[k]);
                }
            }
        }
        delete[] roots;
              
        response[i] = pointRespM4(i, tRoots[0], alpha, templates);
        orientation[i] = tRoots[0];
       
        double temp;
        for (int k=1;k<nr;k++) {
            temp = pointRespM4(i, tRoots[k], alpha, templates);
            if (temp > response[i]) {
                response[i] = temp;
                orientation[i] = tRoots[k];
            }
        }
        delete[] z;
        delete[] tRoots;
    }
}



void filterM5(double** templates, int nx, int ny, double* alpha, double* response, double* orientation) {
    
    double* gx = templates[0];
    double* gy = templates[1];
    double* gxxx = templates[2];
    double* gxxy = templates[3];
    double* gxyy = templates[4];
    double* gyyy = templates[5];
    double* gxxxxx = templates[6];
    double* gxxxxy = templates[7];
    double* gxxxyy = templates[8];
    double* gxxyyy = templates[9];
    double* gxyyyy = templates[10];
    double* gyyyyy = templates[11];

    /*double a11 = alpha[0];
    double a31 = alpha[1];
    double a33 = alpha[2];
    double a51 = alpha[3];
    double a53 = alpha[4];*/
    
    double a10 = alpha[0];
    double a30 = alpha[1];
    double a32 = alpha[2];
    double a52 = alpha[3];
    double a54 = alpha[4];
        
    /*double a1 = 2.0*a31 - 3.0*a33;
    double a2 = 4.0*a51 - 3.0*a53;
    double a3 = 6.0*a33 - 7.0*a31;
    double a4 = 12.0*a51 - 17.0*a53;
    double a5 = 6.0*a53 - 13.0*a51;
    double a6 = 3.0*a33 - 5.0*a31;
    double a7 = a31 - 3.0*a33;
    double a8 = 30.0*a53 - 34.0*a51;
    double a53x2 = 2.0*a53;
    double a11x2 = 2.0*a11;*/
    
    double A, B, C, D, E, F;
    int nr, nt, deg;
    double delta;
    
    for (int i=0;i<nx*ny;++i) {
        
        /*A = -a11*gy[i] + a1*gxxy[i] - a31*gyyy[i] + a53x2*gxxxxy[i] + a2*gxxyyy[i] - a51*gyyyyy[i];
        B = -a11*gx[i] + a1*gxxx[i] + a53x2*gxxxxx[i] + a3*gxyy[i] + a4*gxxxyy[i] + a5*gxyyyy[i];
        C = -a11x2*gy[i] + a6*gxxy[i] + a4*gxxxxy[i] + a7*gyyy[i] + a8*gxxyyy[i] + a2*gyyyyy[i];
        D = -a11x2*gx[i] + a7*gxxx[i] + a2*gxxxxx[i] + a6*gxyy[i] + a8*gxxxyy[i] + a4*gxyyyy[i];
        E = -a11*gy[i] + a3*gxxy[i] + a5*gxxxxy[i] + a1*gyyy[i] + a4*gxxyyy[i] + a53x2*gyyyyy[i];
        F = -a11*gx[i] + a1*gxyy[i] - a31*gxxx[i] - a51*gxxxxx[i]  + a2*gxxxyy[i] + a53x2*gxyyyy[i];*/
                
        /*  sinT5 * (-a10*gx[i] + (2.0*a32-3.0*a30)*gxyy[i] - a32*gxxx[i] + (4.0*a54-3.0*a52)*gxxxyy[i] + 2.0*a52*gxyyyy[i] -a54*gxxxxx[i])
        + cosT5 * ( a10*gy[i] + (3.0*a30-2.0*a32)*gxxy[i] + a32*gyyy[i] + (3.0*a52-4.0*a54)*gxxyyy[i] - 2.0*a52*gxxxxy[i]+ a54*gyyyyy[i])
        + cosT*sinT4 * ( a10*gy[i] + (3.0*a30-2.0*a32)*gyyy[i] + (7.0*a32-6.0*a30)*gxxy[i] - 2.0*a52*gyyyyy[i] + (17.0*a52-12.0*a54)*gxxyyy[i] + (13.0*a54-6.0*a52)*gxxxxy[i])
        + cosT4*sinT * (-a10*gx[i] + (2.0*a32-3.0*a30)*gxxx[i] + (6.0*a30-7.0*a32)*gxyy[i] + 2.0*a52*gxxxxx[i] + (12.0*a54-17.0*a52)*gxxxyy[i] + (6.0*a52-13.0*a54)*gxyyyy[i])       
        + cosT3*sinT2 * ( 2.0*a10*gy[i] + (5.0*a32-3.0*a30)*gxxy[i] + (3.0*a30-a32)*gyyy[i] + (17.0*a52-12.0*a54)*gxxxxy[i] + (34.0*a54-30.0*a52)*gxxyyy[i] + (3.0*a52-4.0*a54)*gyyyyy[i])
        + cosT2*sinT3 * (-2.0*a10*gx[i] + (3.0*a30-5.0*a32)*gxyy[i] + (a32-3.0*a30)*gxxx[i] + (12.0*a54-17.0*a52)*gxyyyy[i] + (30.0*a52-34.0*a54)*gxxxyy[i] + (4.0*a54-3.0*a52)*gxxxxx[i]);*/
        
        A = -a10*gx[i] + (2.0*a32-3.0*a30)*gxyy[i] - a32*gxxx[i] + (4.0*a54-3.0*a52)*gxxxyy[i] + 2.0*a52*gxyyyy[i] -a54*gxxxxx[i];
        B = a10*gy[i] + (3.0*a30-2.0*a32)*gyyy[i] + (7.0*a32-6.0*a30)*gxxy[i] - 2.0*a52*gyyyyy[i] + (17.0*a52-12.0*a54)*gxxyyy[i] + (13.0*a54-6.0*a52)*gxxxxy[i];
        C = -2.0*a10*gx[i] + (3.0*a30-5.0*a32)*gxyy[i] + (a32-3.0*a30)*gxxx[i] + (12.0*a54-17.0*a52)*gxyyyy[i] + (30.0*a52-34.0*a54)*gxxxyy[i] + (4.0*a54-3.0*a52)*gxxxxx[i];
        D = 2.0*a10*gy[i] + (5.0*a32-3.0*a30)*gxxy[i] + (3.0*a30-a32)*gyyy[i] + (17.0*a52-12.0*a54)*gxxxxy[i] + (34.0*a54-30.0*a52)*gxxyyy[i] + (3.0*a52-4.0*a54)*gyyyyy[i];
        E = -a10*gx[i] + (2.0*a32-3.0*a30)*gxxx[i] + (6.0*a30-7.0*a32)*gxyy[i] + 2.0*a52*gxxxxx[i] + (12.0*a54-17.0*a52)*gxxxyy[i] + (6.0*a52-13.0*a54)*gxyyyy[i];
        F = a10*gy[i] + (3.0*a30-2.0*a32)*gxxy[i] + a32*gyyy[i] + (3.0*a52-4.0*a54)*gxxyyy[i] - 2.0*a52*gxxxxy[i]+ a54*gyyyyy[i];
     
        
        A = approxZero(A);
        B = approxZero(B);
        C = approxZero(C);
        D = approxZero(D);
        E = approxZero(E);
        F = approxZero(F);
        
        double* roots;
        double* z;
        
        if (A == 0.0) { // quartic
            if (B == 0.0) { // cubic
                if (C == 0.0) { // quadratic
                    if (D == 0.0) { // linear
                        if (E == 0.0) { // null
                            deg = 1;
                            z = new double[2*deg];
                            z[0] = 0.0;
                            z[1] = 0.0;
                        } else {
                            deg = 1;
                            double a[2] = {E, D};
                            z = new double[2*deg];
                            gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                            gsl_poly_complex_solve(a, deg+1, w, z);
                            gsl_poly_complex_workspace_free(w);
                        }
                    } else { // solve quadratic
                        deg = 2;
                        double a[3] = {E, D, C};
                        z = new double[2*deg];
                        gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                        gsl_poly_complex_solve(a, deg+1, w, z);
                        gsl_poly_complex_workspace_free(w);
                    }
                } else { // solve cubic
                    if ( (D == 0.0) && (F == 0.0) ) {
                        delta = -E/C;
                        if (delta > 0.0) {
                            deg = 3;
                            delta = sqrt(delta);
                            z = new double[2*deg];
                            z[0] = 0.0;
                            z[1] = 0.0;
                            z[2] = delta;
                            z[3] = 0.0;
                            z[4] = -delta;
                            z[5] = 0.0;
                        } else {
                            deg = 1;
                            z = new double[2*deg];
                            z[0] = 0.0;
                            z[1] = 0.0;
                        }
                    } else {
                        deg = 3;
                        double a[4] = {F, E, D, C};
                        z = new double[2*deg];
                        gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                        gsl_poly_complex_solve(a, deg+1, w, z);
                        gsl_poly_complex_workspace_free(w);
                    }
                }
            } else { // solve quartic
                deg = 4;
                double a[5] = {F, E, D, C, B};
                z = new double[2*deg];
                gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
                gsl_poly_complex_solve(a, deg+1, w, z);
                gsl_poly_complex_workspace_free(w);
            }
        } else {
            deg = 5;
            double a[6] = {F, E, D, C, B, A};
            z = new double[2*deg];
            gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(deg+1);
            gsl_poly_complex_solve(a, deg+1, w, z);
            gsl_poly_complex_workspace_free(w);
        }

        
        // # real roots
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0)
                nr++;
        }
        roots = new double[nr];
        nr = 0;
        for (int k=0;k<deg;++k) {
            if (z[2*k+1]==0.0) {
                roots[nr++] = z[2*k];
            }
        }

        
        double* tRoots;
        if (nr == 0) { //orientation[i] = 0.0;
            nt = 4;
            tRoots = new double[nt];
            tRoots[0] = -PI/2.0;
            tRoots[1] = 0.0;
            tRoots[2] = PI/2.0;
            tRoots[3] = PI;
        } else {
            nt = 2*nr;
            tRoots = new double[nt];
            for (int k=0;k<nr;k++) {
                tRoots[k] = atan(roots[k]);
                tRoots[k+nr] = opposite(tRoots[k]);
            }
        }
        delete[] roots;
              
        response[i] = pointRespM5(i, tRoots[0], alpha, templates);
        orientation[i] = tRoots[0];
       
        double temp;
        for (int k=1;k<nt;k++) {
            temp = pointRespM5(i, tRoots[k], alpha, templates);
            if (temp > response[i]) {
                response[i] = temp;
                orientation[i] = tRoots[k];
            }
        }
        delete[] z;
        delete[] tRoots;
    }
}



// Mirror position in image domain for interpolation
int mirror(int x, int nx) {
    if (x >= 0 && x < nx) {
        return x;
    } else if (x < 0) {
        return -x;
    } else {
        return 2*nx-2-x;
    }
}


double interp(double* image, int nx, int ny, double x, double y) {
    int xi = (int)x;
    int yi = (int)y;
    int x0, x1, y0, y1;
    
    double dx = x-xi;
    double dy = y-yi;
    if (x < 0) { dx = -dx; x1 = mirror(xi-1, nx); } else { x1 = mirror(xi+1, nx); }
    if (y < 0) { dy = -dy; y1 = mirror(yi-1, ny); } else { y1 = mirror(yi+1, ny); }
    x0 = mirror(xi, nx);
    y0 = mirror(yi, ny);
    return (1.0-dy)*(dx*image[x1+y0*nx] + (1.0-dx)*image[x0+y0*nx]) + dy*(dx*image[x1+y1*nx] + (1.0-dx)*image[x0+y1*nx]);
}



void computeNMS(double* response, double* orientation, double* nms, int nx, int ny) {
    
    double ux, uy, v1, v2;
    
    div_t divRes;
    for (int i=0;i<nx*ny;++i) {
        divRes = div(i, nx);
        ux = cos(orientation[i]);
        uy = sin(orientation[i]);        
        v1 = interp(response, nx, ny, divRes.rem+ux, divRes.quot+uy);
        v2 = interp(response, nx, ny, divRes.rem-ux, divRes.quot-uy);
        if (v1 > response[i] || v2 > response[i]) {
            nms[i] = 0.0;
        } else {
            nms[i] = response[i];
        }        
    }
}






void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        
    // Input checks
    if (nrhs != 3)
        mexErrMsgTxt("Required inputs arguments: image, filter order, sigma.");
    
    if (nlhs > 4)
        mexErrMsgTxt("Too many output arguments");
    
    if (!mxIsDouble(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2)
        mexErrMsgTxt("Input image must be a 2-D array.");
    
    if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1 || *mxGetPr(prhs[1])<1 || *mxGetPr(prhs[1])>5)
        mexErrMsgTxt("The order 'M' must be an integer between 1 and 5.");
    int M = (int) *mxGetPr(prhs[1]);
    
    if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 || *mxGetPr(prhs[2]) <= 0.0)
        mexErrMsgTxt("Sigma must be a strictly positive scalar value.");
    double sigma = *mxGetPr(prhs[2]);
        
    size_t nx = mxGetN(prhs[0]); // cols
    size_t ny = mxGetM(prhs[0]);
    int N = nx*ny;
    
    int L = 2*(int)(4.0*sigma)+1; // support of the Gaussian kernels
    
    if (L>nx || L>ny) {
        mexPrintf("Sigma must be smaller than %.2f\n", (fmin(nx,ny)-1)/8.0);
        mexErrMsgTxt("Sigma value results in filter support that is larger than image.");
    }
    
    double* input = mxGetPr(prhs[0]);
    double* pixels = (double*)malloc(sizeof(double)*N);
    
    // Process inputs
    
    // Switch matrix to row-major (Matlab uses column-major)
    div_t divRes;
    for (int i=0;i<N;++i) {
        divRes = div(i, ny);
        pixels[divRes.quot+divRes.rem*nx] = input[i];
    }
    
    
    // number of partial derivative templates
    int nTemplates = getTemplateN(M);
    
    
    // Allocate template memory
    double** templates = new double*[nTemplates];
    for (int i=0;i<nTemplates;++i) {
        templates[i] = new double[N];
    }
    double* response = new double[N];
    double* orientation = new double[N];
    
    // Compute the templates used in the steerable filterbank
    computeBaseTemplates(pixels, nx, ny, M, sigma, templates);
    double* alpha = getWeights(M, sigma);

    
    // apply filter
    switch (M) {
        case 1:
            filterM1(templates, nx, ny, alpha, response, orientation);
            break;
        case 2:
            filterM2(templates, nx, ny, alpha, response, orientation);
            break;
        case 3:
            filterM3(templates, nx, ny, alpha, response, orientation);
            break;
        case 4:
            filterM4(templates, nx, ny, alpha, response, orientation);
            break;
        case 5:
            filterM5(templates, nx, ny, alpha, response, orientation);
            break;
    }
   
    
    // Process outputs
    // Switch outputs back to column-major format
    
    if (nlhs > 0) {
        for (int i=0;i<N;++i) {
            divRes = div(i, nx);
            pixels[divRes.quot+divRes.rem*ny] = response[i];
        }
        plhs[0] = mxCreateDoubleMatrix(ny, nx, mxREAL);
        memcpy(mxGetPr(plhs[0]), pixels, N*sizeof(double));
    }
    
    if (nlhs > 1) { // return orientation map
        for (int i=0;i<N;++i) {
            divRes = div(i, nx);
            pixels[divRes.quot+divRes.rem*ny] = orientation[i];
        }
        plhs[1] = mxCreateDoubleMatrix(ny, nx, mxREAL);
        memcpy(mxGetPr(plhs[1]), pixels, N*sizeof(double));
    }
    
    if (nlhs > 2) { // Apply NMS
        double* nms = new double[N];
        computeNMS(response, orientation, nms, nx, ny);
        for (int i=0;i<N;++i) {
            divRes = div(i, nx);
            pixels[divRes.quot+divRes.rem*ny] = nms[i];
        }
        delete[] nms;
        plhs[2] = mxCreateDoubleMatrix(ny, nx, mxREAL);
        memcpy(mxGetPr(plhs[2]), pixels, N*sizeof(double));
    }
    
    if (nlhs > 3) { // return filterbank
        const mwSize dims[] = {ny, nx, nTemplates};
        plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* p = mxGetPr(plhs[3]);
        
        for (int z=0;z<nTemplates;++z) {
            for (int i=0;i<N;++i) {
                divRes = div(i, nx);
                p[z*N + divRes.quot+divRes.rem*ny] = templates[z][i];
            }
        }
        /*int nx = 100;
        int ny = nx;
        double sigma = 4.0;
        int w = ceil(4.0*sigma);        
        plhs[3] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        double* gx = mxGetPr(plhs[3]);
        double theta = PI/3;
        for (int x=0;x<nx;++x) {
            for (int y=0;y<ny;++y) {
                gx[y+x*ny] = -(x-w)/sigma/sigma*exp(-((x-w)*(x-w)+(y-w)*(y-w))/(2.0*sigma*sigma))*cos(theta)
                    -(y-w)/sigma/sigma*exp(-((x-w)*(x-w)+(y-w)*(y-w))/(2.0*sigma*sigma)) * sin(theta);
            }
        */
        
        
        
        
    }    
    
    
    // Free memory
    for (int i=0;i<nTemplates;++i) {
        delete [] templates[i];
    }
    delete[] templates;
    delete[] response;
    delete[] orientation;
    delete[] alpha;
}


// compile with:
// export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2010b.app/bin/maci64 && g++ -ansi -Wall -g -DARRAY_ACCESS_INLINING -I. -L/Applications/MATLAB_R2010b.app/bin/maci64 -I../../mex/include/ -I/Applications/MATLAB_R2010b.app/extern/include steerableDetector.cpp -lmx -lmex -lgsl
// test with:
// valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./a.out 2>&1 | grep steerable

/*int main(void) {
    int nx = 200;
    int ny = 100;
    int N = nx*ny;
    double* pixels = new double[nx*ny];
    for (int i=0;i<nx*ny;++i) {
        pixels[i] = rand();
    }
    int M = 3;
    double sigma = 3.0;
    
    
    int nTemplates = getTemplateN(M);

    
    // Allocate template memory
    double** templates = new double*[nTemplates];
    for (int i=0;i<nTemplates;++i) {
        templates[i] = new double[N];
    }
    double* response = new double[N];
    double* orientation = new double[N];
    
    // Compute the templates used in the steerable filterbank
    computeBaseTemplates(pixels, nx, ny, M, sigma, templates);
    double* alpha = getWeights(M, sigma);

    
    // apply filter
    switch (M) {
        case 1:
            filterM1(templates, nx, ny, alpha, response, orientation);
            break;
        case 2:
            filterM2(templates, nx, ny, alpha, response, orientation);
            break;
        case 3:
            filterM3(templates, nx, ny, alpha, response, orientation);
            break;
        case 4:
            filterM4(templates, nx, ny, alpha, response, orientation);
            break;
        case 5:
            filterM5(templates, nx, ny, alpha, response, orientation);
            break;
    }
    
    for (int i=0;i<nTemplates;++i) {
        delete [] templates[i];
    }
    delete[] templates;
    delete[] response;
    delete[] orientation;
    delete[] alpha;
    
    delete[] pixels;
}*/

