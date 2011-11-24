#include <stdio.h>
#include <math.h>
#include "riskregression.h"


void comptestfunc(times,Ntimes,px,cu,vcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antpers,gamma,line,timepow)
double *times,*cu,*vcu,*vcudif,*test,*testOBS,*Ut,*simUt,*gamma,*timepow;
int *px,*Ntimes,*antsim,*weighted,*antpers,*line;
matrix **W4t;
{
  matrix *Delta,*tmpM1;
  vector *gammavt,*tmpv1t,*tmpv1,*rowX,*xi,*difX,*ssrow,*VdB,
    *gammai[*antpers],*gammav; 
  /*float gasdev(),expdev(),ran1(); 
   */
  double norm_rand(); 
  void GetRNGstate(),PutRNGstate();  
  int i,k,l,s,c;
  double xij,vardif,tau,time,dtime,random,fabs(),sqrt(),stime,mtime;// unused var:x
  double *cumweight=calloc(*px,sizeof(double));

  malloc_vecs(*px,&tmpv1t,&tmpv1,&rowX,&xi,&difX,&ssrow,&VdB,&gammavt,&gammav,NULL); 
  malloc_mat(*Ntimes,*px,Delta); malloc_mat(*Ntimes,*px,tmpM1);
  for (i=0;i<*antpers;i++) malloc_vec(*px,gammai[i]); 

  /* printf("Simulations start N= %ld \n",*antsim);  */

  GetRNGstate();  /* to use R random normals */

  stime=times[0];
  mtime=times[(*Ntimes-1)]-stime; 
  tau=times[(*Ntimes-1)]-times[0]; Ut[0]=times[0]; 
  if (*weighted==0) {
    if (*line==0) for (i=0;i<*px;i++) cumweight[i]=tau;
    else for (i=0;i<*px;i++) cumweight[i]=mtime*mtime*0.5-stime*stime*0.5;}

  /* computation of constant effects */ 
  for (i=0;i<*px;i++) {
    if (fabs(timepow[i])<0.01) { // timepow ca 0
      for (s=0;s<*Ntimes;s++)
      {  // time=times[s];dtime=times[s]-times[s-1];
	 if (vcu[(i+1)*(*Ntimes)+s]>0) {
	 cumweight[i]=cumweight[i]+(1/vcu[(i+1)*(*Ntimes)+s]);
	 gamma[i]=gamma[i]+cu[(i+1)*(*Ntimes)+s]/vcu[(i+1)*(*Ntimes)+s];

	 for (c=0;c<*antpers;c++) VE(gammai[c],i)=
	     VE(gammai[c],i)+ME(W4t[c],s,i)/vcu[(i+1)*(*Ntimes)+s]; 
	 }
      } 
      gamma[i]=gamma[i]/cumweight[i]; 
      VE(gammav,i)=gamma[i]; 
      for (c=0;c<*antpers;c++) VE(gammai[c],i)=VE(gammai[c],i)/cumweight[i]; 
    } else  { 
      gamma[i]=cu[(i+1)*(*Ntimes)+(*Ntimes-1)]/pow(mtime,timepow[i]);; 
      VE(gammav,i)=gamma[i]; 
      for (c=0;c<*antpers;c++) VE(gammai[c],i)=ME(W4t[c],*Ntimes-1,i)/pow(mtime,timepow[i]);; 
    }
  } /*  i=1..px */


  /*
    if (*weighted>=1) {
    for (s=1;s<*Ntimes;s++) {
    vec_zeros(VdB);

    for (i=0;i<*antpers;i++) {
    extract_row(W4t[i],s,tmpv1);  
    vec_subtr(tmpv1,gammai[i],difX);
    vec_star(difX,difX,rowX);
    vec_add(rowX,VdB,VdB); }
    for (k=1;k<=*px;k++) vcudif[k*(*Ntimes)+s]=VE(VdB,k-1); 
    }
    } 
  */ /* weighted==1 */ 

  scl_vec_mult(1,gammav,gammavt); 
  /* Computation of observed teststatistics */ 
  for (s=1;s<*Ntimes;s++)
    {
 if (vcu[i*(*Ntimes)+s]>0) {
      time=times[s]-times[0];dtime=times[s]-times[s-1];
   
      for (i=1;i<=*px;i++) {
	xij=fabs(cu[i*(*Ntimes)+s])/sqrt(vcu[i*(*Ntimes)+s]);

	if (xij>testOBS[i-1]) testOBS[i-1]=xij; } 

      for (i=1;i<=*px;i++) VE(xi,i-1)=cu[i*(*Ntimes)+s];
//      if (*line==1) scl_vec_mult(time,gammav,gammavt); 
      for (i=0;i<*px;i++) VE(gammavt,i)=VE(gammav,i)*pow(time,timepow[i]);; 
      vec_subtr(xi,gammavt,difX); vec_star(difX,difX,ssrow); 

      Ut[s]=times[s]; 

      for (i=0;i<*px;i++) { 
        if (*weighted>=2) vardif=vcudif[(i+1)*(*Ntimes)+s];  else vardif=1; 
	if (*weighted>=2)  {
	  if ((s>*weighted) && (s<*Ntimes-*weighted))  
	    VE(difX,i)=VE(difX,i)/sqrt(vardif); else VE(difX,i)=0.0; 
	} else VE(difX,i)=VE(difX,i); 

        Ut[(i+1)*(*Ntimes)+s]=VE(difX,i);
        c=(*px)+i;
	if (fabs(VE(difX,i))>testOBS[c]) testOBS[c]=fabs(VE(difX,i));
        c=2*(*px)+i;
        if ((s>*weighted) && (s<*Ntimes-*weighted))  
	  testOBS[c]=testOBS[c]+VE(ssrow,i)*dtime; }
    } 
    }
  /*  for (i=0;i<3*(*px);i++) printf(" %lf \n",testOBS[i]);  */

  /* simulation of testprocesses and teststatistics */ 
  for (k=1;k<*antsim;k++) { 
    mat_zeros(Delta); vec_zeros(tmpv1); 
    for (i=0;i<*antpers;i++) {
      /* random=gasdev(&idum);  */
      random=norm_rand();
      scl_mat_mult(random,W4t[i],tmpM1);mat_add(tmpM1,Delta,Delta); 
      scl_vec_mult(random,gammai[i],xi); vec_add(xi,tmpv1,tmpv1);}
      scl_vec_mult(1,tmpv1,tmpv1t); 

    for (s=1;s<*Ntimes;s++) { 
    if (vcu[i*(*Ntimes)+s]>0) {
      time=times[s]-times[0]; dtime=times[s]-times[s-1]; 
      extract_row(Delta,s,rowX); 
//      if (*line==1) scl_vec_mult(times[s],tmpv1,tmpv1t); 
      for (i=0;i<*px;i++) VE(tmpv1t,i)=VE(tmpv1,i)*pow(time,timepow[i]);; 
      vec_subtr(rowX,tmpv1t,difX); vec_star(difX,difX,ssrow); 

      for (i=0;i<*px;i++) { 
	VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(vcu[(i+1)*(*Ntimes)+s]);
	if (VE(xi,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(xi,i); 
	if (*weighted>=1) vardif=vcudif[(i+1)*(*Ntimes)+s];  else vardif=1; 
	if (*weighted>=2)  {
	  if ((s>*weighted) && (s<*Ntimes-*weighted))  
	    VE(difX,i)=VE(difX,i)/sqrt(vardif); else VE(difX,i)=0; 
	} else VE(difX,i)=VE(difX,i); 

        if (k<51) {l=(k-1)*(*px)+i; simUt[l*(*Ntimes)+s]=VE(difX,i);}

        c=(*px+i); VE(difX,i)=fabs(VE(difX,i));
        if (VE(difX,i)>test[c*(*antsim)+k]) test[c*(*antsim)+k]=VE(difX,i);
	c=2*(*px)+i; 
        if ((s>*weighted) && (s<*Ntimes-*weighted))  
	  test[c*(*antsim)+k]=test[c*(*antsim)+k]+VE(ssrow,i)*dtime/vardif; 
      }
    } }  /* s=1..Ntimes */ 
  }  /* k=1..antsim */ 

  PutRNGstate();   /* to use R random normals */

  free_mats(&Delta,&tmpM1,NULL); 
  free_vecs(&gammavt,&tmpv1t,&VdB,&rowX,&difX,&xi,&tmpv1,&ssrow,&gammav,NULL); 
  for (i=0;i<*antpers;i++) free_vec(gammai[i]);  
  free(cumweight); 
}
