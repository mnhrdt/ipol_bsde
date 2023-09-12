#include <cmath>
#include <ctime>
#include <omp.h>
class dbCart{
        public:
        float x;
        float y;
        dbCart(float,float);
};
class dbGrad{
	public:
		float n;
		float dx;
		float dy;
};
class dbRGB{
	public:
		float R;
		float G;
		float B;
		dbRGB(float,float,float);
		dbRGB(float);
		dbRGB();		
		friend dbRGB operator+(const dbRGB &, const dbRGB &);
		friend dbRGB operator-(const dbRGB &, const dbRGB &);
		friend dbRGB operator*(const dbRGB &, const dbRGB &);
		friend dbRGB operator/(const dbRGB &, const dbRGB &);	
		static float abs(dbRGB);
		static dbRGB round(dbRGB);
		static dbRGB max(float,dbRGB);
		static dbRGB min(float,dbRGB);
};
class NormalDist{
        private:
                        float* N;
                        const static long int MAX = 10000000;
                        int position;
        public:
                        NormalDist();
                        float Next();
			~NormalDist();
};
    class dbImage
	{
		public:
			float* A;
			int Width;
			int Height;
			size_t getWidth();
			size_t getHeight();
			float* getArray();
			dbImage();
			dbImage(size_t,size_t);
			dbImage(float*,size_t,size_t);
			float& operator()(int,int);
			float operator()(float,float);
			~dbImage();
			dbGrad grad(float,float,float);
};
class dbImageRGB
	{
		public:
			dbRGB* A;
			int Width;
			int Height;
			size_t getWidth();
			size_t getHeight();
			dbRGB* getArray();
			dbImageRGB();
			dbImageRGB(size_t,size_t);
			dbImageRGB(dbRGB*,size_t,size_t);
			dbRGB& operator()(int,int);
			dbRGB operator()(float,float);
			~dbImageRGB();
			dbGrad grad(float,float,float);		
};
double PSNR(dbImageRGB&, dbImageRGB&);
double PSNR(dbImage&, dbImage&);
