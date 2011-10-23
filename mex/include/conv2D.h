/* 2-D convolution using folded (symmetric/anti-symmetric) kernels
 *
 * (c) Francois Aguet, 2008
 * Adapted from Java (ijtools package) Jul 12, 2011
 * */


void convolveEvenX(double *input, int nx, int ny, double *kernel, int k, double* output) {
    
    int k_1 = k-1;
    int b = 2*nx-2;
    
    int idx = 0;
    for (int y=0;y<ny;y++) {
        for (int x=0;x<k_1;x++) {
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<=x;i++) {
                output[idx] += kernel[i]*(input[idx-i]+input[idx+i]);
            }
            for (int i=x+1;i<k;i++) {
                output[idx] += kernel[i]*(input[i-x+y*nx]+input[idx+i]);
            }
            idx++;
        }
        for (int x=k_1;x<=nx-k;x++) {
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<k;i++) {
                output[idx] += kernel[i]*(input[idx-i]+input[idx+i]);
            }
            idx++;
        }
        for (int x=nx-k_1;x<nx;x++) {
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<nx-x;i++) {
                output[idx] += kernel[i]*(input[idx-i]+input[idx+i]);
            }
            for (int i=nx-x;i<k;i++) {
                output[idx] += kernel[i]*(input[b-i-x+y*nx]+input[idx-i]);
            }
            idx++;
        }
    }
}


void convolveEvenY(double *input, int nx, int ny, double *kernel, int k, double* output) {
    
    int k_1 = k-1;
    int b = 2*ny-2;
    
    int idx, inx;
    for (int x=0;x<nx;x++) {
        for (int y=0;y<k_1;y++) {
            idx = x+y*nx;
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<=y;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]+input[idx+inx]);
            }
            for (int i=y+1;i<k;i++) {
                output[idx] += kernel[i]*(input[(i-y)*nx+x]+input[idx+i*nx]);
            }
        }
        for (int y=k_1;y<=ny-k;y++) {
            idx = x+y*nx;
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<k;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]+input[idx+inx]);
            }
        }
        for (int y=ny-k_1;y<ny;y++) {
            idx = x+y*nx;
            output[idx] = kernel[0]*input[idx];
            for (int i=1;i<ny-y;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]+input[idx+inx]);
            }
            for (int i=ny-y;i<k;i++) {
                output[idx] += kernel[i]*(input[(b-i-y)*nx+x]+input[idx-i*nx]);
            }
        }
    }
}


void convolveOddX(double *input, int nx, int ny, double *kernel, int k, double* output) {
    
    int k_1 = k-1;
    int b = 2*nx-2;
    
    int idx = 0;
    for (int y=0;y<ny;y++) {
        for (int x=0;x<k_1;x++) {
            output[idx] = 0.0;
            for (int i=1;i<=x;i++) {
                output[idx] += kernel[i]*(input[idx-i]-input[idx+i]);
            }
            for (int i=x+1;i<k;i++) {
                output[idx] += kernel[i]*(input[i-x+y*nx]-input[idx+i]);
            }
            idx++;
        }
        for (int x=k_1;x<=nx-k;x++) {
            output[idx] = 0.0;
            for (int i=1;i<k;i++) { // value at 0 is 0
                output[idx] += kernel[i]*(input[idx-i]-input[idx+i]); // conv -> flip
            }
            idx++;
        }
        for (int x=nx-k_1;x<nx;x++) {
            output[idx] = 0.0;
            for (int i=1;i<nx-x;i++) {
                output[idx] += kernel[i]*(input[idx-i]-input[idx+i]);
            }
            for (int i=nx-x;i<k;i++) {
                output[idx] += kernel[i]*(input[idx-i]-input[b-i-x+y*nx]);
            }
            idx++;
        }
    }
}


void convolveOddY(double *input, int nx, int ny, double *kernel, int k, double* output) {
    
    int k_1 = k-1;
    int b = 2*ny-2;
    
    int idx, inx;
    for (int x=0;x<nx;x++) {
        for (int y=0;y<k_1;y++) {
            idx = x+y*nx;
            output[idx] = 0.0;
            for (int i=1;i<=y;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]-input[idx+inx]);
            }
            for (int i=y+1;i<k;i++) {
                output[idx] += kernel[i]*(input[(i-y)*nx+x]-input[idx+i*nx]);
            }
        }
        for (int y=k_1;y<=ny-k;y++) {
            idx = x+y*nx;
            output[idx] = 0.0;
            for (int i=1;i<k;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]-input[idx+inx]);
            }
        }
        for (int y=ny-k_1;y<ny;y++) {
            idx = x+y*nx;
            output[idx] = 0.0;
            for (int i=1;i<ny-y;i++) {
                inx = i*nx;
                output[idx] += kernel[i]*(input[idx-inx]-input[idx+inx]);
            }
            for (int i=ny-y;i<k;i++) {
                output[idx] += kernel[i]*(input[idx-i*nx]-input[(b-i-y)*nx+x]);
            }
        }
    }
}
