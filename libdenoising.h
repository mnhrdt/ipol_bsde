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
#include <iostream>
#include "io_png.h"
#include "lib.h"

/*------Auxiliary function displaying the progress of the algorithm---------*/
void printProgress (double);

/*---Declaration of the function c(t) for greyscale and RGB images----------*/ 
float function_c(float,float);

/*---Declaration of the function b(t) for greyscale and RGB images----------*/ 
float function_b(float);

/*---Declaration of the function g() for greyscale images-------------------*/ 
float function_g(int,int,dbImage&,int,int,float,float,float);

/*---Declaration of the function g() for RGB images-------------------------*/ 
float function_g(int,int,dbImageRGB&,int,int,float,float,float);

/*---Declaration of the main function for RGB denoising---------------------*/ 
/* Arguments: initial image, the noise, noisy image, denoised image,--------*/ 
/*------------add noise option, enhancing parameter-------------------------*/
/* Return: the PSNR---------------------------------------------------------*/
double BSDEdenoisingRGB(dbImageRGB&,float,dbImageRGB&,dbImageRGB&,int,float);

/*---Declaration of the main function for greyscale  denoising--------------*/ 
/* Arguments: initial image, the noise, noisy image, denoised image,--------*/ 
/*------------add noise option, enhancing parameter-------------------------*/
/* Return: the PSNR---------------------------------------------------------*/
double BSDEdenoising(dbImage&,float,dbImage&,dbImage&,int,float);
