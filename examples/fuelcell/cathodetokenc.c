#include<stdio.h>
#include<string.h>
#include<stdlib.h>

static int init=1;
static double R,T,m1,m2,m3,kappa,eta,Pref,epsilon;
static double D12tilde,D13tilde,D23tilde,D21tilde,D31tilde,D32tilde;
static double x1c,x2c,x3c,pc,RT,m13,m23;

// ######################################################################

int equal(char *s, int length, char *t)
{

  int i,len1,len2;

  for (i=0;i<length;i++){
    if (s[i]==' ') break;
  }

  len1=i;

  for (i=0;i<1000;i++){
    if (t[i]=='\0') break;
  }

  len2=i;

  if (len1==len2){
    for (i=0;i<len1;i++)
      if (s[i]!=t[i]) return 0;
      
    return 1;
  }
  else
    return 0;
}

// ######################################################################

int getline(char *s,int limit,FILE *fp)
{
  int i;
  char c;

  i=0;
  while (i<limit){
    c=fgetc(fp);
    if (c==EOF || c=='\n') break;
    else{
      s[i]=c;
      i++;
    }
  }

  s[i]='\0';

  if (c==EOF)
    return -1;
  else
    return i;   
}

// ######################################################################

int split(char *s, char *left, char *right)
{
  int eq,i,j,k;

  eq=-1;

  i=0; j=0; k=0;

  while (1){
    if (s[i]=='\0') break;
    
    if (s[i]=='=') eq=i;
    else if (s[i]!=' ')
      if (eq==-1){
	left[j]=s[i];
	j++;
      }
      else {
	right[k]=s[i];
	k++;
      }
    i++;
  }

  left[j]='\0';
  right[k]='\0';

  return eq;
}

// ######################################################################

void NLinit()
{
  // ----------------------------------------------------------------------
  // If the first call, read from cathode.con
  // ----------------------------------------------------------------------
  
  int i;
  char line[1000],lefts[1000],rights[1000];
  FILE *fp;

  if (init==1){
		 
    printf("This is the C version. Zinc linked to cathodetokenc.c\n");
    printf("NLinit, reading cathodetokenc.con\n");

    init=0;

    fp=fopen("cathodetoken.con","r");
  
    while (getline(line,1000,fp)!=-1){

      i=split(line,lefts,rights);

      if (i!=-1) {
           
	if (strcmp(lefts,"R")==0)
	  sscanf(rights,"%lf", &R);

	else if (strcmp(lefts,"T")==0)
	  sscanf(rights,"%lf",&T);
             
	else if (strcmp(lefts,"m1")==0)
	  sscanf(rights,"%lf",&m1);
              
	else if (strcmp(lefts,"m2")==0)
	  sscanf(rights,"%lf",&m2);
              
	else if (strcmp(lefts,"m3")==0)
	  sscanf (rights,"%lf",&m3);
              
	else if (strcmp(lefts,"kappa")==0)
	  sscanf (rights,"%lf",&kappa);
              
	else if (strcmp(lefts,"eta")==0)
	  sscanf (rights,"%lf",&eta);
              
	else if (strcmp(lefts,"Pref")==0)
	  sscanf (rights,"%lf",&Pref);
              
	else if (strcmp(lefts,"epsilon")==0)
	  sscanf (rights,"%lf",&epsilon);
              
	else if (strcmp(lefts,"D12tilde")==0)
	  sscanf (rights,"%lf",&D12tilde);
	else if (strcmp(lefts,"D13tilde")==0)
	  sscanf (rights,"%lf",&D13tilde);
	else if (strcmp(lefts,"D23tilde")==0)
	  sscanf (rights,"%lf",&D23tilde);
              
	else if (strcmp(lefts,"x1c")==0)
	  sscanf (rights,"%lf",&x1c);
	else if (strcmp(lefts,"x2c")==0)
	  sscanf (rights,"%lf",&x2c);
	else if (strcmp(lefts,"x3c")==0)
	  sscanf (rights,"%lf",&x3c);
	else if (strcmp(lefts,"pc")==0)
	  sscanf (rights,"%lf",&pc);

      }
    }

    printf("pc=%f\n",pc);

    fclose(fp);

    m13=m1-m3;
    m23=m2-m3;
    RT=R*T;

    D21tilde=D12tilde;
    D31tilde=D13tilde;
    D32tilde=D23tilde;
  }
}

// ######################################################################

double cfun(char *label, double *x, double *y, double *z, 
       double ur[3], double dur[3][3], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
				       
  double D11,D12,D13,D21,D22,D23,D31,D32,D33;
  double x1,x2,p,x3,w1,w2,w3,m,denom;

  // ----------------------------------------------------------------------
  // Set values appropriately
  // ----------------------------------------------------------------------

  NLinit();
  
  x1=ur[0];
  x2=ur[1];
  p=ur[2];

  x3=1-x1-x2;

  m=m1*x1+m2*x2+m3*x3;

  w1=m1*x1/m;
  w2=m2*x2/m;
  w3=m3*x3/m;

  denom=x1/(D12tilde*D13tilde)+x2/(D12tilde*D23tilde)+x3/(D13tilde*D23tilde);

  D11=((w2+w3)*(w2+w3)/(x1*D23tilde)+w2*w2/(x2*D13tilde)+w3*w3/(x3*D12tilde))/denom;
  D22=((w1+w3)*(w1+w3)/(x2*D13tilde)+w1*w1/(x1*D23tilde)+w3*w3/(x3*D12tilde))/denom;
  D33=((w1+w2)*(w1+w2)/(x3*D12tilde)+w1*w1/(x1*D23tilde)+w2*w2/(x2*D13tilde))/denom;

  D12=(w1*(w2+w3)/(x1*D23tilde)+w2*(w1+w3)/(x2*D13tilde)-w3*w3/(x3*D12tilde))/denom;
  D13=(w1*(w2+w3)/(x1*D23tilde)+w3*(w1+w2)/(x3*D12tilde)-w2*w2/(x2*D13tilde))/denom;
  D23=(w2*(w1+w3)/(x2*D13tilde)+w3*(w1+w2)/(x3*D12tilde)-w1*w1/(x1*D23tilde))/denom;

  D21=D12;
  D31=D13;
  D32=D23;

  if (equal(label,length,"$B11"))
    return p*m1*x1*(D11-D13)/RT;

  else if (equal(label,length,"$B12"))
    return p*m1*x1*(D12-D13)/RT;

  else if (equal(label,length,"$B21"))
    return p*m2*x2*(D21-D23)/RT;

  else if (equal(label,length,"$B22"))
    return p*m2*x2*(D22-D23)/RT;

  else if (equal(label,length,"$Q1"))
    return m1*x1/RT*(x1*(D11*(1-m1/m)-D13*(1-m3/m)) 
		   + x2*(D12*(1-m2/m)-D13*(1-m3/m))  
	                             +D13*(1-m3/m));

  else if (equal(label,length,"$Q2"))
    return m2*x2/RT*(x1*(D21*(1-m1/m)-D23*(1-m3/m))
                   + x2*(D22*(1-m2/m)-D23*(1-m3/m)) 
                                     +D23*(1-m3/m));

  else if (equal(label,length,"$E"))
    return p*m*kappa/(RT*eta);
  else {
    printf("Error in Cfun: token not recognised\n");
    exit(1);
  }
  
}

// ######################################################################

double ffun(char *label, double *x, double *y, double *z, 
       double ur[3], double dur[3][3], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{

  double x1,x2,p,x3,m,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz;
  double mx,my,mz;

  NLinit();

  x1=ur[0];
  x2=ur[1];
  p=ur[2];

  x3=1-x1-x2;

  m=m1*x1+m2*x2+m3*x3;

  x1x=dur[0][0];
  x1y=dur[1][0];
  x1z=dur[2][0];

  x2x=dur[0][1];
  x2y=dur[1][1];
  x2z=dur[2][1];
  
  px=dur[0][2];
  py=dur[1][2];
  pz=dur[2][2];

  mx=m1*x1x + m2*x2x + m3*(-x1x-x2x);
  my=m1*x1y + m2*x2y + m3*(-x1y-x2y);
  mz=m1*x1z + m2*x2z + m3*(-x1z-x2z);

  if (equal(label,length,"$F1"))
    return p*kappa*m1/(eta*RT)*(px*(x1x-x1*mx/m)+py*(x1y-x1*my/m)+pz*(x1z-x1*mz/m));
  else if (equal(label,length,"$F2"))
    return p*kappa*m2/(eta*RT)*(px*(x2x-x2*mx/m)+py*(x2y-x2*my/m)+pz*(x2z-x2*mz/m));
  else {
    printf ("ffun: unknown token encountered\n");
    exit(1);
  }

}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void seturp(double urp[3], double ur[3])
{
  urp[0]=ur[0];
  urp[1]=ur[1];
  urp[2]=ur[2];
}

// ######################################################################

void setdurp(double durp[3][3],double dur[3][3])
{
  int i,j;

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      durp[i][j]=dur[i][j];
}

// ######################################################################

void dcfun_du(char *token, double *x, double *y, double *z,
	      double ur[3], double dur[3][3], int *nvar,int *istep,
	      int *ireg, int *iregup, double *rnode, double *vec,
	      int *imax, int *jmax, int *kmax, double dCdu[3], int length)
{

  double urp[3];
  double x1cent,x2cent,pcent,dx1,dx2,dp;

  x1cent =0.5;    // x1,x2=[0,1]
  x2cent =0.5;
  pcent =111430;   // near 1 atm

  dx1=x1cent*1e-3;
  dx2=x2cent*1e-3;
  dp=  pcent*1e-3;

  seturp(urp,ur);
  urp[0]=urp[0]+dx1;
  dCdu[0] = (cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	    -cfun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dx1;

  seturp(urp,ur);
  urp[1]=urp[1]+dx2;
  dCdu[1] = (cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
            -cfun(token,x,y,z,ur ,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dx2;

  seturp(urp,ur);
  urp[2]=urp[2]+dp;
  dCdu[2] = (cfun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	    -cfun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dp;
}

// ######################################################################
     

void dcfun_ddu(char *token, double *x, double *y, double *z,
	       double ur[3], double dur[3][3], int *nvar,int *istep,
	       int *ireg, int *iregup, double *rnode, double *vec,
	       int *imax, int *jmax, int *kmax, double dCddu[3][3], int length)
{
  double durp[3][3];
  double x1cent,x1range,x2cent,x2range,pcent,prange,dist,ddx1,ddx2,ddp;
  int i,j;

  x1cent =0.5;    // x1,x2=[0,1]
  x1range=0.1;

  x2cent =0.5;
  x2range=0.1;

  pcent =111430;   // near 1 atm
  prange= 10000;

  dist=2e-3;

  ddx1=x1range/dist*0.0001;
  ddx2=x2range/dist*0.0001;
  ddp=  prange/dist*0.0001;

  // x1
  for(j=0;j<3;j++){
    setdurp(durp,dur);
    durp[j][0] += ddx1;
    dCddu[j][0] = (cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
		  -cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddx1;
  }

  // x2
  for(j=0;j<3;j++){
    setdurp(durp,dur);
    durp[j][1] += ddx2;
    dCddu[j][1] = (cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
		  -cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddx2;
  }

  // p
  for(j=0;j<3;j++){
    setdurp(durp,dur);
    durp[j][2] += ddp;
    dCddu[j][2] = (cfun(token,x,y,z,ur,durp,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
         	  -cfun(token,x,y,z,ur,dur,   nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddp;
  }
}

// ######################################################################

void dffun_du(char *token, double *x, double *y, double *z,
	      double ur[3], double dur[3][3], int *nvar,int *istep,
	      int *ireg, int *iregup, double *rnode, double *vec,
	      int *imax, int *jmax, int *kmax, double dfdu[3], int length)
{

  double urp[3];
  double x1cent,x1range,x2cent,x2range,pcent,prange,dist,dx1,dx2,dp;

  double x1,x2,p,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz,topx,topy,topz,bot;
  static int analytical=1;

  NLinit();

  if (analytical){
    x1=ur[0];
    x2=ur[1];
    p=ur[2];

    x1x=dur[0][0];  //(1,1);
    x1y=dur[1][0];  //(1,2);
    x1z=dur[2][0];  //(1,3);

    x2x=dur[0][1];  //(2,1);
    x2y=dur[1][1];  //(2,2);
    x2z=dur[2][1];  //(2,3);
  
    px=dur[0][2];   //(3,1);
    py=dur[1][2];   //(3,2);
    pz=dur[2][2];   //(3,3);

    topx=(m1-m3)*x1x+(m2-m3)*x2x;
    topy=(m1-m3)*x1y+(m2-m3)*x2y;
    topz=(m1-m3)*x1z+(m2-m3)*x2z;
     
    bot =(m1-m3)*x1 +(m2-m3)*x2 +m3;

    if (equal(token,length,"$F1")) {
        dfdu[0]=p*kappa*m1/(eta*RT)
             *(px*(-topx/bot+x1*topx/bot*bot*(m1-m3))
             + py*(-topy/bot+x1*topy/bot*bot*(m1-m3))
	     + pz*(-topz/bot+x1*topz/bot*bot*(m1-m3)));

        dfdu[1]=p*kappa*m1/(eta*RT)
             *(px*(          x1*topx/bot*bot*(m2-m3))
             + py*(          x1*topy/bot*bot*(m2-m3))
	     + pz*(          x1*topz/bot*bot*(m2-m3)));

        dfdu[2]=kappa*m1/(eta*RT)
             *(px*(x1x-x1*topx/bot)
             + py*(x1y-x1*topy/bot) 
	     + pz*(x1z-x1*topz/bot));
    }
    else if (equal(token,length,"$F2")){

        dfdu[0]=p*kappa*m2/(eta*RT)
             *(px*(          x2*topx/bot*bot*(m1-m3))
             + py*(          x2*topy/bot*bot*(m1-m3))
             + pz*(          x2*topz/bot*bot*(m1-m3)));

        dfdu[1]=p*kappa*m2/(eta*RT)
             *(px*(-topx/bot+x2*topx/bot*bot*(m2-m3))
             + py*(-topy/bot+x2*topy/bot*bot*(m2-m3))
	       + pz*(-topz/bot+x2*topz/bot*bot*(m2-m3)));

        dfdu[2]=kappa*m2/(eta*RT)
             *(px*(x2x-x2*topx/bot) 
             + py*(x2y-x2*topy/bot) 
	       + pz*(x2z-x2*topz/bot));
    }
    else {
      printf("dffun_du: Unknown token\n");
      exit(1);
    }
  }
  else {

    x1cent =0.5;    // x1,x2=[0,1]
    x2cent =0.5;
    pcent =111430;   // near 1 atm
     
    dx1=x1cent*1e-3;
    dx2=x2cent*1e-3;
    dp=  pcent*1e-3;

    seturp(urp,ur);
    urp[0] += dx1;
    dfdu[0]=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	    -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dx1;
     
    seturp(urp,ur);
    urp[1] += dx2;
    dfdu[1]=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	    -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dx2;

    seturp(urp,ur);
    urp[2] += dp;
    dfdu[2]=(ffun(token,x,y,z,urp,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	    -ffun(token,x,y,z,ur, dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/dp;     
  }
}

// ######################################################################

void dffun_ddu(char *token, double *x, double *y, double *z,
	       double ur[3], double dur[3][3], int *nvar,int *istep,
	       int *ireg, int *iregup, double *rnode, double *vec,
	       int *imax, int *jmax, int *kmax, double dfddu[3][3], int length)
{
  double durp[3][3];
  double x1cent,x1range,x2cent,x2range,pcent,prange,dist,ddx1,ddx2,ddp;
  int i,j;

  double x1,x2,p,x1x,x1y,x1z,x2x,x2y,x2z,px,py,pz,bot;
  static int analytical=1;

  NLinit();

  if (analytical){
    x1=ur[0];
    x2=ur[1];
    p=ur[2];

    x1x=dur[0][0];  //(1,1);
    x1y=dur[1][0];  //(1,2);
    x1z=dur[2][0];  //(1,3);

    x2x=dur[0][1];  //(2,1);
    x2y=dur[1][1];  //(2,2);
    x2z=dur[2][1];  //(2,3);
  
    px=dur[0][2];   //(3,1);
    py=dur[1][2];   //(3,2);
    pz=dur[2][2];   //(3,3);

    bot=(m1-m3)*x1+(m2-m3)*x2+m3;

    if (equal(token,length,"$F1")) {
      dfddu[0][0]=p*kappa*m1/(eta*RT)*(1-x1*(m1-m3)/bot)*px;
      dfddu[1][0]=p*kappa*m1/(eta*RT)*(1-x1*(m1-m3)/bot)*py;
      dfddu[2][0]=p*kappa*m1/(eta*RT)*(1-x1*(m1-m3)/bot)*pz;

      dfddu[0][1]=p*kappa*m1/(eta*RT)*( -x1*(m2-m3)/bot)*px;
      dfddu[1][1]=p*kappa*m1/(eta*RT)*( -x1*(m2-m3)/bot)*py;
      dfddu[2][1]=p*kappa*m1/(eta*RT)*( -x1*(m2-m3)/bot)*pz;

      dfddu[0][2]=p*kappa*m1/(eta*RT)*(x1x-x1*((m1-m3)*x1x+(m2-m3)*x2x)/bot);
      dfddu[1][2]=p*kappa*m1/(eta*RT)*(x1y-x1*((m1-m3)*x1y+(m2-m3)*x2y)/bot);
      dfddu[2][2]=p*kappa*m1/(eta*RT)*(x1z-x1*((m1-m3)*x1z+(m2-m3)*x2z)/bot);
    }
    else if (equal(token,length,"$F2")) {
      dfddu[0][0]=p*kappa*m2/(eta*RT)*( -x2*(m1-m3)/bot)*px;
      dfddu[1][0]=p*kappa*m2/(eta*RT)*( -x2*(m1-m3)/bot)*py;
      dfddu[2][0]=p*kappa*m2/(eta*RT)*( -x2*(m1-m3)/bot)*pz;

      dfddu[0][1]=p*kappa*m2/(eta*RT)*(1-x2*(m2-m3)/bot)*px;
      dfddu[1][1]=p*kappa*m2/(eta*RT)*(1-x2*(m2-m3)/bot)*py;
      dfddu[2][1]=p*kappa*m2/(eta*RT)*(1-x2*(m2-m3)/bot)*pz;

      dfddu[0][2]=p*kappa*m2/(eta*RT)*(x2x-x2*((m1-m3)*x1x+(m2-m3)*x2x)/bot);
      dfddu[1][2]=p*kappa*m2/(eta*RT)*(x2y-x2*((m1-m3)*x1y+(m2-m3)*x2y)/bot);
      dfddu[2][2]=p*kappa*m2/(eta*RT)*(x2z-x2*((m1-m3)*x1z+(m2-m3)*x2z)/bot);
    }
    else {
      printf("dffun_ddu: bad token\n");
      exit(1);
    }
  }
  else {

    x1cent =0.5;    // x1,x2=[0,1]
    x1range=0.1;
     
    x2cent =0.5;
    x2range=0.1;
     
    pcent =111430;   // near 1 atm
    prange= 10000;
     
    dist=2e-3;
     
    ddx1=x1range/dist*0.0001;
    ddx2=x2range/dist*0.0001;
    ddp=  prange/dist*0.0001;
     
    // x1
    for(j=0;j<3;j++){
      setdurp(durp,dur);
      durp[j][0]+=ddx1;
      dfddu[j][0]=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	          -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddx1;
    }

    // x2
    for(j=0;j<3;j++){
      setdurp(durp,dur);
      durp[j][1]+=ddx2;
      dfddu[j][1]=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
		  -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddx2;
    }

     // p
    for(j=1;j<3;j++){
      setdurp(durp,dur);
      durp[j][2]+=ddp;
      dfddu[j][2]=(ffun(token,x,y,z,ur,durp, nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length)
	          -ffun(token,x,y,z,ur,dur,  nvar,istep,ireg,iregup,rnode,vec,imax,jmax,kmax,length))/ddp;
    }
  }

}


