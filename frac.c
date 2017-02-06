#include <FPT.h>
#include <stdio.h>
#include <math.h>

int width=1200,height,maxIter=1000;
double SHSL, L;

typedef struct{
 double re;
 double im;
}COMPLEX;

double mod(COMPLEX z){
 //return sqrt(z.re*z.re + z.im*z.im);
 return z.re*z.re+z.im*z.im;
}


//int setColor(int n, double mag){
// double nsmooth = n+1-log(log(mag))/log(2);
// double C = (1-abs(2*L-1))*SHSL;
// double H = nsmooth/60;
// double X = C*(1-abs(H%2-1));
// 
//
// G_rgb(r,g,b);
// return 1;
//}

int draw(double reMin, double reMax, double imMin, double imMax){
 double shade;
 int i,j,iteration;
 COMPLEX c,z,temp;
 for(i=0;i<width;i++){
  for(j=0;j<height;j++){
    c.re = (((double)i)/width)*abs(reMax-reMin)+reMin;
    c.im = (((double)j)/height)*abs(imMax-imMin)+imMin;
    z.re = 0;
    z.im = 0;
    iteration = 0;
    while((mod(z)<=2*2) && (iteration<maxIter)){
	temp.re = z.re*z.re - z.im*z.im + c.re;
	temp.im = 2*z.re*z.im + c.im;
	
	z.re = temp.re;
	z.im = temp.im;
	iteration = iteration +1;
    }

    if(iteration == maxIter){
	G_rgb(0,0,0);
    }
    else{
	shade = (iteration-log(log(mod(z))/log(2)))/iteration;
	G_rgb(.5,shade,.5);
    }
    G_point(i,j);
  }
 }
}


int main(){
 height = (int) ((2/3.5)*width);
 double p[2],reMin=-2.5,reMax=1,imMin=-1,imMax=1;
 double centerx=0,centery=0,rspan,ispan;
 G_init_graphics(width,height);
 while(1){
  G_rgb(0,0,0);
  G_clear();
  draw(reMin,reMax,imMin,imMax);
  G_wait_click(p);

  rspan = reMax-reMin;
  ispan = imMax-imMin;
  centerx = (p[0]/3.5)*(rspan)+reMin;
  centery = (p[1]/2)*(ispan)+imMin;

  reMin = centerx-(rspan*.3);
  reMax = centerx+(rspan*.3);
  imMin = centery-(ispan*.3);
  imMax = centery+(ispan*.3);
 }
 double q=G_wait_key();
 G_close();
}

