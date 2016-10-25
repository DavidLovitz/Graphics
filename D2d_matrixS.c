#include <D2d_matrix.h>

/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( 1 )          (1)

instead of (x',y',1) = (x,y,1) * M  

*/



int D2d_print_mat (double a[3][3])
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 




int D2d_copy_mat (double a[3][3], double b[3][3])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 




int D2d_make_identity (double a[3][3])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
}


double dot_product(double a[3], double b[3]){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}



int D2d_mat_mult (double res[3][3], double a[3][3], double b[3][3]){
  double new_mat[3][3],temp1[3],temp2[3];
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
	temp1[0]=a[r][0];temp1[1]=a[r][1];temp1[2]=a[r][2];
	temp2[0]=b[0][c];temp2[1]=b[1][c];temp2[2]=b[2][c];
	new_mat[r][c] = dot_product(temp1,temp2);
      }
  }
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
	res[r][c] = new_mat[r][c];
      }
  }
  return 1;
}





/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


int D2d_translate (double a[3][3], double b[3][3], double dx, double dy)
// a = translation*a  
// b = b*translation_inverse  
{
  double t[3][3] ;

  D2d_make_identity(t) ;

  t[0][2] =  dx ;  t[1][2] = dy ;  
  D2d_mat_mult(a,  t,a) ;

  t[0][2] = -dx ;  t[1][2] = -dy ;
  D2d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D2d_scale (double a[3][3], double b[3][3], double sx, double sy)
// a = scale*a  
// b = b*scale_inverse
{
  double t[3][3];
  D2d_make_identity(t);

  //store scaling matrix in a
  t[0][0] = sx; t[1][1] = sy;
  D2d_mat_mult(a, t,a);

  //avoid division by 0
  if(sx==0 || sy==0){
    return 0;
  }
  //store inverse in b
  t[0][0] = 1/sx; t[1][1] = 1/sy;
  D2d_mat_mult(b, t,b);

  return 1;
}


int D2d_rotate (double a[3][3], double b[3][3], double radians)
// a = rotate*a  
// b = b*rotate_inverse
{
  double t[3][3];
  D2d_make_identity(t);

  //store rotation matrix in a
  t[0][0] = cos(radians); t[0][1] = -sin(radians);
  t[1][0] = sin(radians); t[1][1] = cos(radians);
  D2d_mat_mult(a, t,a);

  //store inverse in b
  t[0][0] = cos(-radians); t[0][1] = -sin(-radians);
  t[1][0] = sin(-radians); t[1][1] = cos(-radians);
  D2d_mat_mult(b, t,b);

  return 1;
}


int D2d_negate_x (double a[3][3], double b[3][3])
// negate the x....reflects in the y-axis
// a = reflect*a 
// b = b*reflect_inverse  
{
  double t[3][3];
  D2d_make_identity(t);

  //store scaling matrix in a
  t[0][0] = -1; t[1][1] = 1;
  D2d_mat_mult(a, t,a);

  //store inverse in b
  t[0][0] = -1; t[1][1] = 1;
  D2d_mat_mult(b, t,b);

  return 1;
}


int D2d_negate_y (double a[3][3], double b[3][3])
// negate the y....reflects in the x-axis
// a = reflect*a 
// b = b*reflect_inverse
{
  double t[3][3];
  D2d_make_identity(t);

  //store scaling matrix in a
  t[0][0] = 1; t[1][1] = -1;
  D2d_mat_mult(a, t,a);

  //store inverse in b
  t[0][0] = 1; t[1][1] = -1;
  D2d_mat_mult(b, t,b);
  return 1;
}


int D2d_mat_mult_points (double *X, double *Y,
                         double m[3][3],
                         double *x, double *y, int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like D2d_mat_mult_points (x,y, m, x,y, n) ;
{
  double temp;
  int c;
 
  for (c = 0 ; c < numpoints ; c++ ) {
    temp = x[c]*m[0][0]+y[c]*m[0][1] + m[0][2];
    Y[c] = x[c]*m[1][0]+y[c]*m[1][1] + m[1][2];
    X[c] = temp;
  }
  return 1;
}
