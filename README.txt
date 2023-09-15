#TITLE

Implementation of Image Denoising based on Backward Stochastic
Differential Equations


#VERSION

2.0


#AUTHOR

Dariusz Borkowski <dbor@mat.umk.pl>


#CONTENTS

./BSDE.cpp - This program reads PNG images, runs BSDE
	 method and writes PNG images. Only 8bit (grayscale
	 and RGB) PNG images are handled.

./libdenoising.cpp, ./libdenoising.h - Implementation of BSDE
	method: Algorithm 1, Algorithm 2, Algorithm 3.
	This source code provides an implementation of the BSDE
	denoising algorithm.

./lib.cpp, ./lib.h - Implementation of classes dbImageRGB
	and dbImage representing RGB and greyscale images and
	operations on them.

./io_png.c, ./io_png.h - PNG read/write interface by Nicolas
	Limare <nicolas.limare@cmla.ens-cachan.fr>.

./lena.png, ./dice.png, ./traffic.png, ./computer.png,
	./book.png, book10.png - Test PNG images.



The code is written in C/C++ and should compile on any system
with an g++ compiler.

The libpng header and libraries are required on the system for
compilation and execution.

The implementation uses OPENMP which not supported by old
versions of g++.


Ubuntu 22.04 LTS installation:

sudo apt-get install make
sudo apt-get install g++
sudo apt-get install libpng-dev


# COMPILATION INSTRUCTIONS

Simply use the provided makefile, with the command `make`.


# USAGE

usage: ./BSDE original.png sigma noisy.png denoised.png add_noise c

BSDE takes 6 parameter:
* original.png 		: initial noise free image (PNG image)
* sigma     		: the noise standard deviation
				(non negative real number)
* noisy.png  		: noisy image used by the denoising
				algorithm (PNG image)
* denoised.png	 	: denoised image (PNG image)
* add_noise     	: the add noise option (integer number
				from the set {0,1})
* c			: the enhancing parameter c (non-negative
				real number from the range [0, 1])


# EXAMPLES OF USE

Use the provided makefile, with the command `make tests`

Example 1.
We add noise with a standard deviation sigma=20 to lena.png original
image and reconstruct it without the enhancing effect c=0.0. We save
the noisy and reconstructed image as in.png and out.png, respectively.
The program displays the PSNR value of the original and reconstructed
image.

./BSDE lena.png 20 in.png out.png 1 0.0


Example 2.
We reconstruct the already noisy image lena.png using sigma=10 and enhancing
effect c=0.4.  The input image is copied unchanged to file in.png, and the
reconstructed image is saved to file out.png.  Since there is no ground truth
to compare against, the PSNR value is not computed.

./BSDE lena.png 10 in.png out.png 0 0.4


# COPYRIGHT AND LICENSE INFORMATION

Copyright (c) 2023 Dariusz Borkowski

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
