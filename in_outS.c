#include <FPT.h>

int swidth, sheight ;

int click_and_save (double *x, double *y)
{
  int n ;
  double P[2] ;

  G_rgb(0,1,0.5) ;
  G_fill_rectangle(0,0,swidth,20) ;

  G_rgb(1,0,0) ;
  G_wait_click(P);

  n = 0 ;
  while (P[1] > 20) {
    x[n] = P[0] ;
    y[n] = P[1] ;
    G_circle(x[n],y[n],2) ;
    if (n > 0) { G_line(x[n-1],y[n-1], x[n],y[n]) ;}
    n++ ;
    G_wait_click(P) ;
  }

  return n ;
}



int in_out (double *x, double *y, int n, double P[2])
// return 1 if point P is inside the convex polygon
// else return 0
{
  int i;
  double m,cx=0,cy=0,inout,pval;

  for(i=0;i<n;i++){
    cx = cx + x[i];
    cy = cy + y[i];
  }
  cx = cx/n;
  cy = cy/n;

  for(i=0;i<n;i++){
    if(i==n-1){m=(y[i]-y[0])/(x[i]-x[0]);}
    else{m=(y[i]-y[i+1])/(x[i]-x[i+1]);}
    inout = m*(cx-x[i])+y[i]-cy;
    pval = m*(P[0]-x[i])+y[i]-P[1];
    if(inout>0 && pval<0){
      return 0;
    }
    if(inout<0 && pval>0){
      return 0;
    }
  }
  return 1;
}



int main()
{
  double xp[1000],yp[1000] ;
  int n,q,in ;
  double P[2];

  printf("%lf\n",sqrt(10)*sin(atan2(1,3)-M_PI/4));
  swidth = 700 ; sheight = 700 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  G_rgb(1,0,0) ;
  n = click_and_save(xp,yp) ;
  G_rgb(0,1,0) ;
  G_fill_polygon(xp,yp,n) ;
  
 
  while(1){
    G_wait_click(P) ;
    printf("%d\n",in_out(xp,yp,n,P));
    G_rgb(0,0,1) ;
    in = in_out(xp,yp,n,P) ;
    if(in == 1){
      G_rgb(1,0,0);
    }
    G_fill_circle(P[0],P[1],2) ;
  }
  q = G_wait_key() ;
}
