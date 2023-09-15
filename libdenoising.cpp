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

#include "libdenoising.h"

/*---------Auxiliary function displaying the progress of the algorithm-----*/

void printProgress (double percentage){
    int val = (int) (percentage * 100);
    std::cout<< "\r\t"<<val<<"%";
    fflush (stdout);
}

/*---Definition of the function c(t) for greyscale and RGB images-----------*/

float function_c(float t,float c){
	return c;
}

/*---Definition of the function b(t) for greyscale and RGB images-----------*/

float function_b(float t){
	return 0.05f;
}

/*---Definition of the function g() for greyscale images--------------------*/
/*-g() is a measure of neighborhood similarity at points (X1,X2) (x1,x2)-*/

float function_g(int X_1,int X_2,dbImage& Noisy,int x_1,int x_2,float radius,float h,float Sigma){
	float dist=0.0;
	float dif;
    	for (int xx=-radius; xx<= radius; xx++)
        	for (int yy=-radius; yy<=radius; yy++) {
	    		dif=Noisy(X_1+xx,X_2+yy)-Noisy(x_1+xx,x_2+yy);
            		dist += dif*dif;
    		}
	float fdif=std::fmax((dist-2.0f*float(1.0f*(2*radius+1)*(2*radius+1))*Sigma*Sigma),0.0f)/(h*h);
	if (fdif >=  100.0f) return 0.0f;
	return 1.0f-fdif*0.01f;
}

/*---Definition of the function g() for RGB images-------------------------*/
/*-g() is a measure of neighborhood similarity at points (X1,X2) (x1,x2)---*/

float function_g(int X_1,int X_2,dbImageRGB& Noisy,int x_1,int x_2,float radius,float h,float Sigma){
	float dist=0.0;
	dbRGB dif;
    	for (int xx=-radius; xx<= radius; xx++)
        	for (int yy=-radius; yy<=radius; yy++) {
	    		dif=Noisy(X_1+xx,X_2+yy)-Noisy(x_1+xx,x_2+yy);
            		dist += dif.R*dif.R;
            		dist += dif.G*dif.G;
            		dist += dif.B*dif.B;
    		}
	float fdif=std::fmax((dist-2.0f*float(3.0f*(2*radius+1)*(2*radius+1))*Sigma*Sigma),0.0f)/(h*h);
	if (fdif >= 100.0f) return 0.0f;
	return 1.0f-fdif*0.01f;
}

/*--- Implementation for RGB images-----------------------------------------*/

double BSDEdenoisingRGB(dbImageRGB& Original,float Sigma,dbImageRGB& Input,dbImageRGB& totalOutput,int addNoise,float c){

/*---Setting of approximation parameters according to section: 7 Approximation Parameters--*/
	int N=20;
	int j=10+round(sqrtf(Sigma));
	int m=j+3;
	float dt=4.0f;
	float p=Sigma;

/*---Definition of coefficients a(k) according to the formulas (5) (6) and (7)---*/

	float* a=new float[m];
	float tmp;
	for(int k=0;k<=j-1;k++){
		tmp=function_b(dt*k)*dt;
		for(int s=0;s<=k-1;s++)tmp*=(1-function_b(dt*s)*dt);
		a[k]=tmp;
	}
	tmp=1;
	for(int r=j;r<=m-1;r++)tmp*=(1+function_c(dt*j+dt*(r-j),c)*dt);
	tmp=tmp-function_c(dt*j,c)*dt;
	for(int s=0;s<=j-1;s++)tmp*=(1-function_b(dt*s)*dt);
	a[j]=tmp;
	for(int k=j+1;k<=m-1;k++){
		tmp=function_c(dt*j+dt*(k-j),c)*dt;
		for(int r=j;r<=k-1;r++)tmp*=(1+function_c(dt*j+dt*(r-j),c)*dt);
		for(int s=0;s<=j-1;s++)tmp*=(1-function_b(dt*s)*dt);
		a[k]=-tmp;
	}

/*---Setting of parameters for g() function----------------------------------*/

	float radius=1;
	float h=0.55f*Sigma;
     	if(Sigma > 25.0f && Sigma <= 55.0f) { radius = 2; h = 0.4f*Sigma;}
     	if(Sigma > 55.0f) { radius=3; h=0.35*Sigma;}

/*Creating images: noisy image, convolution, gradient norm, partial derivatives*/
/*---The method .Next() is a generator of the normal distribution N(0,1)-------*/

	size_t nx=Original.getWidth();
	size_t ny=Original.getHeight();
    	dbImage totalWeights(nx,ny);
    	dbImageRGB ConvInput(nx,ny);
    	dbImage DX(nx,ny);
    	dbImage DY(nx,ny);
    	dbImage Norm(nx,ny);
	NormalDist gen;
	dbRGB vector;
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			if(addNoise){
				vector.R=gen.Next();
				vector.G=gen.Next();
				vector.B=gen.Next();
				Input(x,y)=Original(x,y)+Sigma*vector;
			}
			else{
				Input(x,y)=Original(x,y);
			}
			totalOutput(x,y).R=0.0f;
			totalOutput(x,y).G=0.0f;
			totalOutput(x,y).B=0.0f;
			totalWeights(x,y)=0.0f;
	}
	for(int x=0;x<(int)nx;x++)
                for(int y=0;y<(int)ny;y++)
                        ConvInput(x,y) = (4 * Input(x, y) + 2 * Input(x + 1, y) + 2 * Input(x - 1, y) +
                        2 * Input(x, y + 1) + 2 * Input(x, y - 1) + Input(x + 1, y + 1) + Input(x - 1, y - 1) +
                        Input(x + 1, y - 1) + Input(x - 1, y + 1)) / float(16.0f);
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
            		dbGrad grad=ConvInput.grad(x,y,0.4);
            		DX(x,y)=grad.dx;
            		DY(x,y)=grad.dy;
            		Norm(x,y)=grad.n;
        }
#ifdef _OPENMP
#pragma omp parallel
#endif
{

/*---Each thread creates private images and variables------------------------------*/

    	dbImageRGB Output(nx,ny);
	dbImage Weights(nx,ny);
	float* X=new float[m];
        float* Y=new float[m];
	float weight;
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			Output(x,y).R=0.0f;
			Output(x,y).G=0.0f;
			Output(x,y).B=0.0f;
			Weights(x,y)=0.0f;
		}

/*------Algorithm 3 = Algorithm1 + Algorithm 2---------------------------------------*/
/*------For each pixel position (x,y) do restoration---------------------------------*/

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int x=0;x<(int)nx;x++){
		printProgress(((double)x+1.0)/int(nx));
		float rand;
		for(int y=0;y<(int)ny;y++){
			X[0]=(float)x;
                	Y[0]=(float)y;

/*---Algorithm 1 -- implementation of Example 2.--------------------------------------*/

if(c>0 && Norm(x,y)>Sigma){
			for(int n=1;n<=N;n++){
				weight=a[0];
				Weights(x,y)+=weight;
				Output(x,y)=Output(x,y)+weight*Input(x,y);
				for(int k=1;k<=j;k++){
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+dt*gen.Next()));
                                        Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+dt*gen.Next()));

/*-----------------Modified diffusion with random terminal time-----------------------*/

					if(dbRGB::abs(ConvInput(X[k],Y[k])-ConvInput(X[k-1],Y[k-1]))>p)k--;
					else {
						weight=a[k];
						Weights(x,y)+=weight;
						Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
					}
				}
				for(int k=j+1;k<m;k++){
					rand=dt*gen.Next();
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+rand*DX(X[k-1],Y[k-1])));
                            		Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+rand*DY(X[k-1],Y[k-1])));
					weight=a[k];
					Weights(x,y)+=weight;
					Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
				}
			}

}

/*---Algorithm 2 -- implementation of Example 3.----------------------------------------*/

else{
			for(int n=1;n<=N;n++){
				weight=a[0];
				Weights(x,y)+=weight;
				Output(x,y)=Output(x,y)+weight*Input(x,y);
				for(int k=1;k<=j;k++){
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+dt*gen.Next()));
                                        Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+dt*gen.Next()));

/*--------------------Modified diffusion with random terminal time------------------------*/

					if(dbRGB::abs(ConvInput(X[k],Y[k])-ConvInput(X[k-1],Y[k-1]))>p)k--;
					else {
						if(k==j){
							weight=a[k];

/*----------------------------------The formula (8)--------------------------------------*/

							for(int i=j+1;i<m;i++)weight+=a[i];
							Weights(x,y)+=weight;
							Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
						}
						else{
						weight=a[k]*function_g(int(round(X[k])),int(round(Y[k])),Input,x,y,radius,h,Sigma);
    						for (int xx=-radius; xx<= radius; xx++)
        						for (int yy=-radius; yy<=radius; yy++) {
								Weights(x+xx,y+yy)+=weight;
								Output(x+xx,y+yy)=Output(x+xx,y+yy)+weight*Input(X[k]+xx,Y[k]+yy);
							}
						}
					}
				}

			}

}
		}

	}
	delete[] X;
	delete[] Y;
#ifdef _OPENMP
#pragma omp critical
#endif
{

/*----------------Gathering partial data from threads---------------------------------*/

	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			totalOutput(x,y)=totalOutput(x,y)+Output(x,y);
			totalWeights(x,y)=totalWeights(x,y)+Weights(x,y);
		}

}
}

/*----------------Normalization of weights--------------------------------------------*/

	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			totalOutput(x,y)=dbRGB::max(0.0f,dbRGB::min(255.0f,totalOutput(x,y)/(totalWeights(x,y))));
		}
	delete[] a;
return PSNR(Original,totalOutput);
}

/*----------------Implementation for greyscale images-----------------------------------*/

double BSDEdenoising(dbImage& Original,float Sigma,dbImage& Input,dbImage& totalOutput,int addNoise,float c){

/*---Setting of approximation parameters according to section: 7 Approximation Parameters--*/

	int N=20;
	int j=10+round(sqrtf(Sigma));
	int m=j+3;
	float dt=4.0f;
	float p=Sigma;

/*---Definition of coefficients a(k) according to the formulas (5) (6) and (7)---------*/

	float* a=new float[m];
	float tmp;
	for(int k=0;k<=j-1;k++){
		tmp=function_b(dt*k)*dt;
		for(int s=0;s<=k-1;s++)tmp*=(1-function_b(dt*s)*dt);
		a[k]=tmp;
	}
	tmp=1;
	for(int r=j;r<=m-1;r++)tmp*=(1+function_c(dt*j+dt*(r-j),c)*dt);
	tmp=tmp-function_c(dt*j,c)*dt;
	for(int s=0;s<=j-1;s++)tmp*=(1-function_b(dt*s)*dt);
	a[j]=tmp;
	for(int k=j+1;k<=m-1;k++){
		tmp=function_c(dt*j+dt*(k-j),c)*dt;
		for(int r=j;r<=k-1;r++)tmp*=(1+function_c(dt*j+dt*(r-j),c)*dt);
		for(int s=0;s<=j-1;s++)tmp*=(1-function_b(dt*s)*dt);
		a[k]=-tmp;
	}

/*---Setting of parameters for g() function----------------------------------*/

	float radius=1;
	float h=0.4f*Sigma;
     	if(Sigma > 15.0f && Sigma <= 30.0f) { radius = 2; h = 0.4f*Sigma;}
     	if(Sigma > 30.0f && Sigma <= 45.0f) { radius = 3; h = 0.35f*Sigma;}
     	if(Sigma > 45.0f && Sigma <= 75.0f) { radius = 4; h = 0.35f*Sigma;}
     	if(Sigma > 75.0f) { radius=5; h=0.3f*Sigma;}

/*Creating images: noisy image, convolution, gradient norm, partial derivatives*/
/*---The method .Next() is a generator of the normal distribution N(0,1)-------*/

	size_t nx=Original.getWidth();
	size_t ny=Original.getHeight();
    	dbImage totalWeights(nx,ny);
    	dbImage ConvInput(nx,ny);
    	dbImage DX(nx,ny);
    	dbImage DY(nx,ny);
    	dbImage Norm(nx,ny);
	NormalDist gen;
	float vector;
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			if(addNoise){
				vector=gen.Next();
				Input(x,y)=Original(x,y)+Sigma*vector;
			}
			else{
				Input(x,y)=Original(x,y);
			}
			totalOutput(x,y)=0.0f;
			totalWeights(x,y)=0.0f;
	}
	for(int x=0;x<(int)nx;x++)
                for(int y=0;y<(int)ny;y++)
                        ConvInput(x,y) = (4 * Input(x, y) + 2 * Input(x + 1, y) + 2 * Input(x - 1, y) +
                        2 * Input(x, y + 1) + 2 * Input(x, y - 1) + Input(x + 1, y + 1) + Input(x - 1, y - 1) +
                        Input(x + 1, y - 1) + Input(x - 1, y + 1)) / float(16.0f);
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
            		dbGrad grad=ConvInput.grad(x,y,0.4);
            		DX(x,y)=grad.dx;
            		DY(x,y)=grad.dy;
			Norm(x,y)=grad.n;
        }
#ifdef _OPENMP
#pragma omp parallel
#endif
{

/*---Each thread creates private images and variables------------------------------*/

    	dbImage Output(nx,ny);
	dbImage Weights(nx,ny);
	float* X=new float[m];
        float* Y=new float[m];
	float weight;
	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			Output(x,y)=0.0f;
			Weights(x,y)=0.0f;
		}

/*------Algorithm 3 = Algorithm1 + Algorithm 2---------------------------------------*/
/*------For each pixel position (x,y) do restoration---------------------------------*/

#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
	for(int x=0;x<(int)nx;x++){
		printProgress(((double)x+1.0)/int(nx));
		float rand;
		for(int y=0;y<(int)ny;y++){
			X[0]=(float)x;
                	Y[0]=(float)y;

/*---Algorithm 1 -- implementation of Example 2.--------------------------------------*/

if(c>0 && Norm(x,y)>Sigma){
			for(int n=1;n<=N;n++){
				weight=a[0];
				Weights(x,y)+=weight;
				Output(x,y)=Output(x,y)+weight*Input(x,y);
				for(int k=1;k<=j;k++){
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+dt*gen.Next()));
                                        Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+dt*gen.Next()));

/*-----------------Modified diffusion with random terminal time-----------------------*/

					if(std::abs(ConvInput(X[k],Y[k])-ConvInput(X[k-1],Y[k-1]))>p)k--;
					else {
						weight=a[k];
						Weights(x,y)+=weight;
						Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
					}
				}
				for(int k=j+1;k<m;k++){
					rand=dt*gen.Next();
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+rand*DX(X[k-1],Y[k-1])));
                            		Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+rand*DY(X[k-1],Y[k-1])));
					weight=a[k];
					Weights(x,y)+=weight;
					Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
				}
			}
}

/*---Algorithm 2 -- implementation of Example 3.----------------------------------------*/

else{
			for(int n=1;n<=N;n++){
				weight=a[0];
				Weights(x,y)+=weight;
				Output(x,y)=Output(x,y)+weight*Input(x,y);
				for(int k=1;k<=j;k++){
					X[k]=std::fmax(0.0f,std::fmin(float(nx),X[k-1]+dt*gen.Next()));
                                        Y[k]=std::fmax(0.0f,std::fmin(float(ny),Y[k-1]+dt*gen.Next()));

/*--------------------Modified diffusion with random terminal time------------------------*/

					if(std::abs(ConvInput(X[k],Y[k])-ConvInput(X[k-1],Y[k-1]))>p)k--;
					else {
						if(k==j){
							weight=a[k];

/*----------------------------------The formula (8)--------------------------------------*/

							for(int i=j+1;i<m;i++)weight+=a[i];
							Weights(x,y)+=weight;
							Output(x,y)=Output(x,y)+weight*Input(X[k],Y[k]);
						}
						else{
						weight=a[k]*function_g(int(round(X[k])),int(round(Y[k])),Input,x,y,radius,h,Sigma);
    						for (int xx=-radius; xx<= radius; xx++)
        						for (int yy=-radius; yy<=radius; yy++) {
								Weights(x+xx,y+yy)+=weight;
								Output(x+xx,y+yy)=Output(x+xx,y+yy)+weight*Input(X[k]+xx,Y[k]+yy);
							}
						}
					}
				}
			}
}
		}

	}
	delete[] X;
	delete[] Y;
#ifdef _OPENMP
#pragma omp critical
#endif
{

/*----------------Gathering partial data from threads---------------------------------*/

	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			totalOutput(x,y)=totalOutput(x,y)+Output(x,y);
			totalWeights(x,y)=totalWeights(x,y)+Weights(x,y);
		}
}
}

/*----------------Normalization of weights--------------------------------------------*/

	for(int x=0;x<(int)nx;x++)
		for(int y=0;y<(int)ny;y++){
			totalOutput(x,y)=std::fmax(0.0f,std::fmin(255.0f,totalOutput(x,y)/(totalWeights(x,y))));
		}
	delete[] a;
	return PSNR(Original,totalOutput);
}

