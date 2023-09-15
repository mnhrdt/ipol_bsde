/*
Copyright (c) 2023 Dariusz Borkowski

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.
*/
#include "lib.h"

dbRGB::dbRGB(float r,float g,float b):R(r),G(g),B(b){}
dbRGB::dbRGB(float f):R(f),G(f),B(f){}
dbRGB::dbRGB():R(0.0),G(0.0),B(0.0){}
dbRGB operator+(const dbRGB &p1, const dbRGB &p2){return dbRGB(p1.R+p2.R,p1.G+p2.G,p1.B+p2.B);}
dbRGB operator-(const dbRGB &p1, const dbRGB &p2){return dbRGB(p1.R-p2.R,p1.G-p2.G,p1.B-p2.B);}
dbRGB operator*(const dbRGB &p1, const dbRGB &p2){return dbRGB(p1.R*p2.R,p1.G*p2.G,p1.B*p2.B);}
dbRGB operator/(const dbRGB &p1, const dbRGB &p2){return dbRGB(p1.R/(p2.R),p1.G/(p2.G),p1.B/(p2.B));}
float dbRGB::abs(dbRGB p){return sqrtf(p.R*p.R+p.G*p.G+p.B*p.B);}
dbRGB dbRGB::round(dbRGB p){return dbRGB(std::round(p.R),std::round(p.G),std::round(p.B));}
dbRGB dbRGB::max(float p1, dbRGB p2){return dbRGB(std::fmax(p1,p2.R),std::fmax(p1,p2.G),std::fmax(p1,p2.B));}
dbRGB dbRGB::min(float p1, dbRGB p2){return dbRGB(std::fmin(p1,p2.R),std::fmin(p1,p2.G),std::fmin(p1,p2.B));}
NormalDist::NormalDist(){
                                N = new float[MAX];
                                srand (time(NULL));
                                float u, v, x;
                                for (int i = 0; i < MAX; i++){
                                        do{
                                                        u = rand()/(float)RAND_MAX;
                                                        v = ((rand()/(float)RAND_MAX) - 0.5) * 2 * sqrtf(2 / expf(1));
                                                        if (u == 0) u = 0.0001;
                                                        x = v / u;
                                                }
                                        while (!(u * u <= expf((-x * x) / 2.0)));
                                        N[i] = x;

                                }
                                position = rand() % (MAX-1);
                        }
float NormalDist::Next(){
                                position++;
                                if (position >= MAX) position = 0;
                                return N[position];
                        }
NormalDist::~NormalDist(){
						delete[] N;
}
size_t dbImage::getWidth(){return Width;}
size_t dbImage::getHeight(){return Height;}
float* dbImage::getArray(){return A;}
dbImage::dbImage(){A=NULL;Width=0;Height=0;}
dbImage::dbImage(size_t W, size_t H){Width = W; Height = H; A = new float[W*H];}
dbImage::dbImage(float* T,size_t W, size_t H){A=T; Width=W; Height=H;}
float& dbImage::operator()(int i,int j){
				if (i < 0) i = 0;
				if (j < 0) j = 0;
				if (i >= (int)Width) i = Width - 1;
				if (j >= (int)Height) j = Height - 1;
				return A[i+j*Width];
			}
float dbImage::operator()(float r,float s){
				float a, b, c, d, deltax, deltay, v;
				int i, j, k, l;
				if (r < 0) r = 0;
				if (s < 0) s = 0;
				if (r >= Width) r = Width - 1;
				if (s >= Height) s = Height - 1;
				i = (int)r;
				j = (int)s;
				k = (i + 1);
				l = (j + 1);
				if (k >= Width) k = Width - 1;
				if (l >= Height) l = Height - 1;
				deltax = (r - i);
				deltay = (s - j);
				a = A[i+j*Width];
				b = A[k+j*Width];
				c = A[i+l*Width];
				d = A[k+l*Width];
				v = (a * (1 - deltax) + b * deltax) * (1 - deltay) + (c * (1 - deltax) + d * deltax) * deltay;
				return v;
			}

dbImage::~dbImage(){
				delete[] A;
			}
dbGrad dbImage::grad(float r, float s, float h){
				float dx, dy;
				dbGrad w;
				dx = (this->operator()(r + h, s) - this->operator()(r - h, s)) / (2 * h);
				dy = (this->operator()(r, s + h) - this->operator()(r, s - h)) / (2 * h);
				w.n = sqrt(dx * dx + dy * dy);
				if (w.n == 0) { w.dx = 0; w.dy = 0; return w;};
				w.dx = dx / w.n;
				w.dy = dy / w.n;
				return w;
}
size_t dbImageRGB::getWidth(){return Width;}
size_t dbImageRGB::getHeight(){return Height;}
dbRGB* dbImageRGB::getArray(){return A;}
dbImageRGB::dbImageRGB(){A=NULL;Width=0;Height=0;}
dbImageRGB::dbImageRGB(size_t W, size_t H){Width = W; Height = H; A = new dbRGB[W*H];}
dbImageRGB::dbImageRGB(dbRGB* T,size_t W, size_t H){A=T; Width=W; Height=H;}
dbRGB& dbImageRGB::operator()(int i,int j){
				if (i < 0) i = 0;
				if (j < 0) j = 0;
				if (i >= (int)Width) i = Width - 1;
				if (j >= (int)Height) j = Height - 1;
				return A[i+j*Width];
			}
dbRGB dbImageRGB::operator()(float r,float s){
				dbRGB a, b, c, d, deltax, deltay, v;
				int i, j, k, l;
				if (r < 0) r = 0;
				if (s < 0) s = 0;
				if (r >= Width) r = Width - 1;
				if (s >= Height) s = Height - 1;
				i = (int)r;
				j = (int)s;
				k = (i + 1);
				l = (j + 1);
				if (k >= Width) k = Width - 1;
				if (l >= Height) l = Height - 1;
				deltax = (r - i);
				deltay = (s - j);
				a = A[i+j*Width];
				b = A[k+j*Width];
				c = A[i+l*Width];
				d = A[k+l*Width];
				v = (a * (1 - deltax) + b * deltax) * (1 - deltay) + (c * (1 - deltax) + d * deltax) * deltay;
				return v;
			}
dbImageRGB::~dbImageRGB(){
			delete[] A;
			}
dbGrad dbImageRGB::grad(float r, float s, float h){
			float dx,dy,Rx,Ry,Gx,Gy,Bx,By,g11,g12,g22,delta,lambda1,norma;
			dbGrad w;
			Rx = ((*this)(r + h, s).R - (*this)(r - h, s).R) / (2.0 * h);
			Ry = ((*this)(r, s + h).R - (*this)(r, s - h).R) / (2.0 * h);
			Gx = ((*this)(r + h, s).G - (*this)(r - h, s).G) / (2.0 * h);
			Gy = ((*this)(r, s + h).G - (*this)(r, s - h).G) / (2.0 * h);
			Bx = ((*this)(r + h, s).B - (*this)(r - h, s).B) / (2.0 * h);
			By = ((*this)(r, s + h).B - (*this)(r, s - h).B) / (2.0 * h);
			g11=Rx*Rx+Gx*Gx+Bx*Bx;
			g12=Rx*Ry+Gx*Gy+Bx*By;
			g22=Ry*Ry+Gy*Gy+By*By;
			delta=((g11-g22)*(g11-g22)+4*g12*g12);
			lambda1=((g11+g22+sqrtf(delta))/(float)2.0);
			/*lambda2=((g11+g22-sqrtf(delta))/(float)2); */
			dx=(2*g12);
			dy=(g22-g11+sqrtf(delta));
			w.n =sqrtf(lambda1/*-lambda2*/);
            norma=sqrtf(dx*dx+dy*dy);
            dx=dx/(norma+0.0001);
            dy=dy/(norma+0.0001);
            if(dx>1)dx=1;
        	if(dy>1)dy=1;
            w.dx = dx ;
            w.dy = dy ;
            return w;
}
double PSNR(dbImageRGB& o, dbImageRGB& r){
	int m=o.getWidth(),n=o.getHeight();
	double mse=0,div=(m*n*3.0);
	for(int i=0;i<m;i++)for(int j=0;j<n;j++){
			mse+=(o(i,j).R-r(i,j).R)*(o(i,j).R-r(i,j).R);
			mse+=(o(i,j).G-r(i,j).G)*(o(i,j).G-r(i,j).G);
			mse+=(o(i,j).B-r(i,j).B)*(o(i,j).B-r(i,j).B);
	}
	mse/=div;
	return 20*log10(255)-10*log10(mse);
}
double PSNR(dbImage& o, dbImage& r){
	int m=o.getWidth(),n=o.getHeight();
	double mse=0,div=(m*n*1.0);
	for(int i=0;i<m;i++)for(int j=0;j<n;j++){
			mse+=(o(i,j)-r(i,j))*(o(i,j)-r(i,j));
	}
	mse/=div;
	return 20*log10(255)-10*log10(mse);
}
