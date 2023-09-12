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

int main(int argc,char **argv){

/*------------------Loading input data---------------------------------------*/

	size_t nx,ny,nc;
    	float *d_v = NULL;
	d_v = io_png_read_f32(argv[1], &nx, &ny, &nc);
	float Sigma=float(atof(argv[2]));
	int addNoise=atoi(argv[5]);
	float c=float(atof(argv[6]));

/*------------------Denoising RGB images----------------------------------*/

if(Sigma>0 && d_v && int(nc)==3){
    	dbRGB* D_V=new dbRGB[nx*ny];
    	for(int i=0;i<(int)(nx*ny);i++){
		D_V[i].R=d_v[i];
		D_V[i].G=d_v[i+nx*ny];
		D_V[i].B=d_v[i+2*nx*ny];
	}
    	dbImageRGB Original(D_V,nx,ny);
    	dbImageRGB Input(nx,ny);
    	dbImageRGB totalOutput(nx,ny);

/*-----------Main algorithm for RGB images--------------------------------*/

	double psnr=BSDEdenoisingRGB(Original,Sigma,Input,totalOutput,addNoise,c);

/*-----------Saving RGB images--------------------------------------------*/

	for(int i=0;i<(int)(nx*ny);i++){
		d_v[i]=Input.A[i].R;
		d_v[i+nx*ny]=Input.A[i].G;
		d_v[i+2*nx*ny]=Input.A[i].B;
	}
	if(io_png_write_f32(argv[3], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";
	for(int i=0;i<(int)(nx*ny);i++){
		d_v[i]=totalOutput.A[i].R;
		d_v[i+nx*ny]=totalOutput.A[i].G;
		d_v[i+2*nx*ny]=totalOutput.A[i].B;
	}
	if(io_png_write_f32(argv[4], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";

/*-----Displaying the PSNR value if we know the original image------------*/

	if(addNoise)std::cout<<"\rPSNR: "<<psnr<<"\n";
	return 0;
}

/*------------------Denoising greyscale images-----------------------------*/

if(Sigma>0 && d_v && int(nc)==1){
    	float* D_V=new float[nx*ny];
    	for(int i=0;i<(int)(nx*ny);i++){
		D_V[i]=d_v[i];
	}
    	dbImage Original(D_V,nx,ny);
    	dbImage Input(nx,ny);
    	dbImage totalOutput(nx,ny);

/*-----------Main algorithm for greyscale images---------------------------*/

	double psnr=BSDEdenoising(Original,Sigma,Input,totalOutput,addNoise,c);

/*-----------Saving greyscale images---------------------------------------*/

	for(int i=0;i<(int)(nx*ny);i++){
		d_v[i]=Input.A[i];
	}
	if(io_png_write_f32(argv[3], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";
	for(int i=0;i<(int)(nx*ny);i++){
		d_v[i]=totalOutput.A[i];
	}
	if(io_png_write_f32(argv[4], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";

/*-----Displaying the PSNR value if we know the original image------------*/

	if(addNoise)std::cout<<"\rPSNR: "<<psnr<<"\n";
	return 0;
}
if(Sigma<=0 && d_v){
	if(io_png_write_f32(argv[3], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";
	if(io_png_write_f32(argv[4], d_v, nx, ny, nc) != 0)std::cout<<"... failed to save png image\n";
	return 0;
}
std::cout<<"error :: not found or not a correct png image\n";
return -1;
}
