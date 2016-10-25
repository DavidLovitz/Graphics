#include <FPT.h> //http://legacy.lclark.edu/~jeff/   
#include <stdio.h>
#include <math.h>
#include <D2d_matrix.h>

//This program displays object.xy style objects
//rotation,scaling,translation, and clipping using a bounding box are all implemented

//How to compile:
//acom -I. D2d_matrixS.c lab3.c

//needs FPT installed at /FPT

//to run: ./a.out obj.xy obj2.xy obj3.xy //etc.

//a and d rotate left and right
//w and s zoom in and out
//clicking sets the vertexes of a bounding box for clipping, must be a convex polygon
//clipping only affects the first object for now

//q quits


//Arrays for storing ojbect's vertexes, rgb components,
//averages and max distances for center of mass
double x[10][3000],y[10][3000],r[10][2000],g[10][2000],
  b[10][2000],xp[10][200][3000],yp[10][200][3000],avgX[10],avgY[10],dmax[10];
//window dimensions, number of objects read from file
  int swidth=600,sheight=600,numobjects;
//vertexes, connectivity data
  int vertexes[10],polygons[10],i,j,k,con[10][2000][20],psize[10][2000];
//bounding information for clipping
double boundx[20], boundy[20],cx,cy;
int bounds=0;


//simple insertion sort for my_fill_polygon
void sort(double *x,int n){
  int s,i,j;
  double t;
  for(j=0;j<n;j++){
    s=j;
    for(i=j+1;i<n;i++){
      if(x[i]<x[s]){s=i;}
    }
    t=x[j]; x[j]=x[s]; x[s]=t;
  }
  return;
}

//my_fill_polygon draws a polygon with n vertexes given by xp and yp
//(not currently being used)
my_fill_polygon(double *xp,double *yp,int n){
  int i,j,k,intercepts;
  double slope,x1,y1,x2,y2;
  for(i=0;i<600;i++){//loop through all y values
    double x[100];
    intercepts = 0;
    for(j=0;j<n;j++){//loop through all segments
      if(j==n-1){//account for final segment
	x1=xp[j];x2=xp[0];y1=yp[j];y2=yp[0];
      }
      else{//store current valeus in shorter variables
	x1 = xp[j];x2=xp[j+1];y1=yp[j];y2=yp[j+1];
      }
      if(((i<y1) && (i>y2)) || ((i>y1) && (i<y2))){//check if i in the right range
	slope = (y1-y2)/(x1-x2);//calculate slope of segment
	x[intercepts]=x1+(i-y1)/slope;//use pt slope form to find x intercept pt w/ y=i
	intercepts++;
      }
    }
    sort(x,intercepts);//sort intercepted x values
    for(k=0;k<intercepts;k=k+2){
      G_line(x[k],i,x[k+1],i);//draw line segments skipping every other (connect,skip,connect,skip,...)
    }
  }
  return;
}

//Stores the bounding window for clipping
//returns the number of vertexes
//max 20 vertexes
int clicked(){
  int n;
  double p[2];

  G_rgb(0,1,0.5);
  G_fill_rectangle(0,0,swidth,20);

  G_rgb(1,1,1);
  G_wait_click(p);

  n=0;
  while(p[1]>20){
    boundx[n] = p[0];
    boundy[n] = p[1];
    G_circle(boundx[n],boundy[n],5);

 
    if(n>0){
      G_line(boundx[n-1],boundy[n-1],boundx[n],boundy[n]);
    }
    n++;
    G_wait_click(p);
  }
  return n;
}
//clips a polygon with n vertexes given by arrays xpt and ypt against the boundx and boundy  polygon
int clip(double *xpt, double *ypt, int n){
  int i, j, k, in[1000], newn=0, nextj, nexti;
  double m, inout, pval, xn[1000], yn[1000], A, B, C, A1, B1, C1;
  cx = 0; cy = 0;
  for(i=0;i<bounds;i++){
    cx += boundx[i];
    cy += boundy[i];
  }
  cx /= bounds;
  cy /= bounds;
  G_rgb(1,1,1);
  G_circle(cx,cy,6);

  //  printf("bounds = %d\n",bounds) ;

  for(i=0;i<bounds;i++){

    if(i==bounds-1){
      nexti = 0;
    }
    else{
      nexti=i+1;
    }
    A=boundy[nexti]-boundy[i];
    B=-boundx[nexti]+boundx[i];
    C=boundy[i]*(boundx[nexti]-boundx[i])-(boundx[i]*(boundy[nexti]-boundy[i]));
    
    inout = A*cx+B*cy+C; 
    
    for(j=0;j<n;j++){
      pval = A*xpt[j] + B*ypt[j] + C;
      if((inout>0 && pval>0) || (inout<0 && pval<0)){
	in[j] = 1;
      }
      else{
	in[j] = 0;
      }
    }
    

    newn = 0 ;
    for(j=0;j<n;j++){
      if(j == n-1){nextj = 0;}
      else{nextj = j+1;}

      A1=ypt[nextj]-ypt[j];
      B1=-xpt[nextj]+xpt[j];
      C1=ypt[j]*(xpt[nextj]-xpt[j])-(xpt[j]*(ypt[nextj]-ypt[j]));

      //if(B1 == 0) printf("problem\n") ;

      if(in[j] && in[nextj]){// in in [store next]
        xn[newn] = xpt[nextj];
	yn[newn] = ypt[nextj];
	newn++;
      }
      if(!in[j] && in[nextj]){//out in [store intersept and next]
	//	xn[newn] = (((B*C1)/B1)-C)/(A-((B*A1)/B1));
	//	yn[newn] = (-(A1*xn[newn])-C1)/B1;
	//CRAMERS RULE to avoid divisions by zero (in most cases)
	xn[newn] = (B*C1 - C*B1) / (A*B1 - B*A1) ;
        yn[newn] = (A1*C - A*C1) / (A*B1 - B*A1) ;


	newn++;
	xn[newn] = xpt[nextj];
	yn[newn] = ypt[nextj];
	newn++;
      }
      if(in[j] && !in[nextj]){//in out [store intercept]
	//	xn[newn] = (((B*C1)/B1)-C)/(A-((B*A1)/B1));
	//	yn[newn] = (-(A1*xn[newn])-C1)/B1;
	xn[newn] = (B*C1 - C*B1) / (A*B1 - B*A1) ;
        yn[newn] = (A1*C - A*C1) / (A*B1 - B*A1) ;

	newn++;
      }
      if(!in[j] && !in[nextj]){//out out
	continue;
      } 
    }
    for(k=0;k<newn;k++){
      xpt[k] = xn[k];
      ypt[k] = yn[k];
    }
    n=newn;
  }
  //new number of vertexes after clipping
  return newn;
}


//draws the object w/ object number onum
void drawobject(int onum){
  int i,j,k;
  double maxX=0,maxY=0,minX=999999,minY=999999;
  double xx[100],yy[100];
  int v,n;
  onum = (int)onum;
  for(i=0;i<polygons[onum];i++){
    //goes through the number of vertexes for polygon j
    n=psize[onum][i];
    for(j=0;j<n;j++){
      v = con[onum][i][j];
      xx[j] = x[onum][v];
      yy[j] = y[onum][v];
    }
    if(bounds>0){
       n=clip(xx,yy,n);
       printf("n = %d ",n) ;
    }
    G_rgb(r[onum][i],g[onum][i],b[onum][i]);//set color
    //my_fill_polygon(xx,yy,n);
    if (n > 0)  G_fill_polygon(xx,yy,n);
  }


  for(i=0;i<bounds;i++){
    G_rgb(1,1,1);
    G_circle(boundx[i],boundy[i],5);
    if(i!=bounds-1){
      G_line(boundx[i],boundy[i],boundx[i+1],boundy[i+1]);
    }
    else{
      G_line(boundx[i],boundy[i],boundx[0],boundy[0]);
    }
  }
}

//reads in from the .xy file all of the object data
void readobject(FILE *f,int onum){
  int i,k;
  //total vertexes
  fscanf(f,"%d",&vertexes[onum]);
  for(i=0;i<vertexes[onum];i++){
    fscanf(f,"%lf %lf",&x[onum][i],&y[onum][i]);
  }
  //total polygons
  fscanf(f,"%d",&polygons[onum]);
  for(i=0;i<polygons[onum];i++){
    //store the number vertexes for polygon i in psize[onum][i]
    fscanf(f,"%d",&psize[onum][i]);
    //go through the number of vertexes and get the vertexes of the polygon i
    //con[onum][i][k] k should go from 0 to psize[onum][i]
    for(k=0;k<psize[onum][i];k++){
      fscanf(f,"%d",&con[onum][i][k]);
    }
  }
  for(i=0;i<polygons[onum];i++){
    fscanf(f,"%lf %lf %lf",&r[onum][i],&g[onum][i],&b[onum][i]);
  }
}

//gets the average x value from object # onum
double getAvgX(int onum){
  double avgX = 0;
  for(i=0; i<vertexes[onum];i++){
    avgX+=x[onum][i];
  }
  return avgX/vertexes[onum];
}
//gets the average y value from object # onum
double getAvgY(int onum){
  double avgY =0;
  for(i=0; i<vertexes[onum];i++){	
    avgY+=y[onum][i];
  }
  return avgY/vertexes[onum];
}

int main(int argc,char **argv)
{
  int i, on=0, clp = 0;
  double key,temp,m[3][3],minverse[3][3];
  numobjects = argc-1;
  //open and read files, store object data
  for(i=0;i<numobjects;i++){
    FILE *f;
    f=fopen(argv[i+1],"r");
    if(f==NULL){printf("file not found");exit(0);}
    readobject(f,i);
  }


  //find avg x and y
  for(i=0;i<numobjects;i++){
    avgX[i] = getAvgX(i);
    avgY[i] = getAvgY(i);
  }


  //find maximum distance from center of mass to outer vertex
  for(i=0;i<numobjects;i++){
    dmax[i] = 0;
    for(j=0;j<vertexes[i];j++){
      temp=sqrt((x[i][j]-avgX[i])*(x[i][j]-avgX[i]) + (y[i][j]-avgY[i])*(y[i][j]-avgY[i]));
      if(temp>dmax[i]){
	dmax[i]=temp;
      }
    }
  }


  //bring objects to origin and scale by 300/dmax then translate to center of screen
  for(i=0;i<numobjects;i++){
    D2d_make_identity(m); D2d_make_identity(minverse);
    D2d_translate(m,minverse,-avgX[i],-avgY[i]);
    D2d_scale(m,minverse,300/dmax[i],300/dmax[i]);
    D2d_translate(m,minverse,300,300);
    D2d_mat_mult_points(x[i],y[i], m, x[i],y[i],vertexes[i]);
 }

  //initialize graphics
  G_init_graphics(swidth, sheight);
  G_rgb(0,0,0);
  G_clear();

  //draw the first object and wait for mouse clicks to define the clipping window 
  drawobject(0);
  bounds = clicked();
  G_rgb(0,0,0);
  G_clear();
  drawobject(0);

  key = G_wait_key();
  //allows for object manipulation until q is pressed
  while(1){
    D2d_make_identity(m); D2d_make_identity(minverse);
    if(key==97){//'a' rotate left
      D2d_translate(m,minverse,-300,-300);
      D2d_rotate(m,minverse,3*M_PI/180);
      D2d_translate(m,minverse,300,300);
    }
    if(key==100){//'d' rotate right
      D2d_translate(m,minverse,-300,-300);
      D2d_rotate(m,minverse,-3*M_PI/180);
      D2d_translate(m,minverse,300,300);
    }
    if(key==119){//'w' zoom in
      D2d_translate(m,minverse,-300,-300);
      D2d_scale(m,minverse,1.2,1.2);
      D2d_translate(m,minverse,300,300);
    }
    if(key==115){//'s' zoom out
      D2d_translate(m,minverse,-300,-300);
      D2d_scale(m,minverse,0.8,0.8);
      D2d_translate(m,minverse,300,300);
    }
    if(key==120){//'x' flip over x axis
      D2d_translate(m,minverse,-300,0);
      D2d_negate_x(m,minverse);
      D2d_translate(m,minverse,300,0);
    }
    if(key==121){//'y' flip over y axis
      D2d_translate(m,minverse,0,-300);
      D2d_negate_y(m,minverse);
      D2d_translate(m,minverse,0,300);
    }
    D2d_mat_mult_points(x[on],y[on], m, x[on],y[on],vertexes[on]);
    if(key<58){on = (int)(key-48);}
    if(key==113){exit(0);}
   
    drawobject(on);
    key=G_wait_key();
    G_rgb(0,0,0);
    G_clear();
  }
  G_close(); // terminate graphics
}
