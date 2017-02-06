#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <FPT.h>

int main(){
 int n=200;
 double theta = 2*M_PI/n, r=300;
 int i,j,k=2;
 G_init_graphics(600,600);
 G_rgb(0,0,0);
 G_clear();
 while(1){
 G_rgb(1,1,1);
 for(i=1;i<n;i++){
   j = (k*i)%n;
   G_rgb(((double)i)/500,.7,.7);
   G_line(r*cos(i*theta)+300,r*sin(i*theta)+300,
	  r*cos(j*theta)+300,r*sin(j*theta)+300);
 }
 double q = G_wait_key();
 k++;
 printf("%d\n",k);
 G_rgb(0,0,0);
 G_clear();
 }
exit(0);
}
