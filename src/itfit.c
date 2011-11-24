#include <stdio.h>
#include <math.h>
#include "riskregression.h"
// {{{ ifit header	 
void itfit(times,
	   Ntimes,
	   x,
	   delta,
	   cause,
	   KMc,
	   z,
	   n,
	   px,
	   Nit,
	   betaS,
	   score,
	   hess,
	   est,
	   var,
	   sim,
	   antsim,
	   test,
	   testOBS,
	   Ut,
	   simUt,
	   weighted,
	   gamma,
	   vargamma,
	   semi,
	   zsem,
	   pg,
	   trans,
	   gamma2,
	   CA,
	   line,
	   detail,
	   biid,
	   gamiid,
	   resample,
	   timepow,
	   clusters,
	   antclust,
	   timepowtest,
	   silent,
	   convc)
     double *times,*betaS,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
     *Ut,*simUt,*gamma,*zsem,*gamma2,*biid,*gamiid,*vargamma,*timepow,
     *timepowtest,*convc;
     int *n,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*weighted,
  *semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,*silent;
{
  // }}}
  // {{{ declare some stuff
  matrix *X,*cX,*A,*AI,*cumAt[*antclust],*VAR,*Z;
  vector *VdB,*risk,*SCORE,*W,*Y,*Gc,*DELTA,*CAUSE,*bhat,*pbhat,*beta,*xi,
    *rr,*rowX,*difbeta,*qs,*bhatub,*betaub,*dcovs,*pcovs,*zi,*rowZ,*zgam; 
  vector *cumhatA[*antclust],*cumA[*antclust],*bet1,*gam,*dp,*dp1,*dp2; 
  int ps,sing,sc,c,i,j,k,l,s,it,convproblems=0; 
  double time,sumscore,totrisk,zgamt, 
    *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double));
  float gasdev(),expdev(),ran1();
  ps=(*px); 
  // }}}
  if (*semi==0) { 
  // {{{ 
    malloc_mat(*n,*px,X); malloc_mat(*n,*px,cX); 
    if (*trans==2) {malloc_mat(*n,*pg,Z);malloc_vecs(*pg,&zgam,&gam,&zi,&rowZ,NULL);}
    malloc_mats(ps,ps,&A,&AI,&VAR,NULL); 

    malloc_vecs(*n,&rr,&bhatub,&risk,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat,NULL); 
    malloc_vecs(*px,&bet1,&xi,&rowX,NULL); 
    malloc_vecs(ps,&dp,&dp1,&dp2,&dcovs,&pcovs,&betaub,&VdB,&qs,&SCORE,&beta,
	       &difbeta,NULL); 

    for (i=0;i<*antclust;i++) {
      malloc_vec(ps,cumhatA[i]); malloc_vec(ps,cumA[i]); 
      malloc_mat(*Ntimes,ps,cumAt[i]);}

    for (c=0;c<ps;c++) VE(beta,c)=betaS[c]; 
    for (c=0;c<*px;c++) VE(bet1,c)=betaS[c]; 
//    if (*trans==0) {for (c=0;c<*pg;c++) VE(gam,c)=betaS[*px+c];}

    for (c=0;c<*n;c++) {VE(Gc,c)=KMc[c]; VE(DELTA,c)=delta[c]; 
      VE(CAUSE,c)=cause[c]; 
      for(j=0;j<*px;j++)  ME(X,c,j)=z[j*(*n)+c]; 
//      if (*trans==0) for(j=0;j<*pg;j++)  ME(Z,c,j)=zsem[j*(*n)+c];
    }

	 
 sc=0;
for (s=0;s<*Ntimes;s++)
{
   time=times[s]; est[s]=time; score[s]=time; var[s]=time;

  for (it=0;it<*Nit;it++)
  {
    totrisk=0; 
    for (j=0;j<*n;j++) { 
      VE(risk,j)=(x[j]>=time); totrisk=totrisk+VE(risk,j);
      extract_row(X,j,xi); 

      VE(bhat,j)=vec_prod(xi,bet1); 

      if (*trans==1) {VE(pbhat,j)=1-exp(-VE(bhat,j));
	scl_vec_mult(1-VE(pbhat,j),xi,dp);}
      if (*trans==2) {
	VE(pbhat,j)=1-exp(-exp(VE(bhat,j))); 
	scl_vec_mult((1-VE(pbhat,j))*exp(VE(bhat,j)),xi,dp); }
      if (*trans==3) {
	    VE(pbhat,j)=exp(VE(bhat,j))/(1+exp(VE(bhat,j))); 
	scl_vec_mult(exp(VE(bhat,j))/pow((1+exp(VE(bhat,j))),2),xi,dp);
      // VE(pbhat,j)=exp(-VE(bhat,j)); scl_vec_mult(-exp(-VE(pbhat,j)),xi,dp);
      }
      if (*trans==4) {
         VE(pbhat,j)=exp(VE(bhat,j)); 
	 scl_vec_mult(exp(VE(bhat,j)),xi,dp);
      }
      if (*trans==8) { // not implemented 
	extract_row(Z,j,zi); vec_star(zi,gam,zgam); zgamt=vec_sum(zgam); 
	VE(rr,j)=exp(zgamt);  
	VE(pbhat,j)=1-(VE(bhat,j)*VE(rr,j))/(1+(VE(bhat,j)*VE(rr,j))); 
	scl_vec_mult((1-VE(pbhat,j))*VE(rr,j),xi,xi); 
	scl_vec_mult((1-VE(pbhat,j))*VE(rr,j)*VE(bhat,j),zi,zi); 
	for (c=0;c<*px;c++) VE(dp,c)=VE(xi,c); 
	for (c=*px;c<ps;c++) VE(dp,c)=VE(zi,c-*px); }
      replace_row(cX,j,dp); 

      VE(Y,j)=((x[j]<=time) & (cause[j]==*CA))*1;
      if (it==*Nit-1) {
	if (KMc[j]<0.00001) vec_zeros(dp); else scl_vec_mult(1/KMc[j],dp,dp); 
	scl_vec_mult(VE(Y,j),dp,dp); vec_add(dp,qs,qs); }
      if (KMc[j]<0.001) VE(Y,j)=(VE(Y,j)/0.001)-VE(pbhat,j); 
      else VE(Y,j)=(VE(Y,j)/KMc[j])-VE(pbhat,j);
    }
    totrisk=vec_sum(risk); MtA(cX,cX,A); 
    invertS(A,AI,silent[0]); sing=0; 
    // head_matrix(cX); print_mat(A); print_mat(AI); 

    if (fabs(ME(AI,0,0))<.0000001) {
      convproblems=1; 
      for (c=0;c<ps;c++) VE(beta,c)=0; 
      for (c=0;c<*px;c++) VE(bet1,c)=betaS[c]; 
      sing=1;
      if (*silent==0) printf("Non-invertible design at time %lf\n",time); 
      it=*Nit-1;  
    }
    if (sing==0) {
      /* print_vec(Y); print_vec(SCORE); print_vec(difbeta); */ 
      vM(cX,Y,SCORE); 
      Mv(AI,SCORE,difbeta); vec_add(beta,difbeta,beta); 
      for (i=0;i<*px;i++) VE(bet1,i)=VE(beta,i); 

      sumscore=0; 
      for (k=0;k<*px;k++) sumscore=sumscore+fabs(VE(difbeta,k)); 
      if ((sumscore<*convc) & (it<*Nit-2)) it=*Nit-2;

      if (isnan(vec_sum(SCORE))) {
	printf("missing values in SCORE %ld \n",(long int) s); 
	convproblems=1; 
	it=*Nit-1; 
	for (c=0;c<ps;c++) VE(beta,c)=0; 
	for (c=0;c<*px;c++) VE(bet1,c)=betaS[c]; 
	}
    }

    if (*detail==1) { 
      printf(" s er %ld, Estimate beta \n",(long int) s); print_vec(beta); 
      printf("Score D l\n"); print_vec(difbeta); 
      printf("Information -D^2 l\n"); print_mat(AI); };

    if (it==*Nit-1) scl_vec_mult(1/totrisk,qs,qs); 
  } /* it */

vec_zeros(VdB); mat_zeros(VAR); 

   for (j=0;j<*antclust;j++) {vec_zeros(cumA[j]);vec_zeros(cumhatA[j]);}
   for (i=0;i<*n;i++) { 
      j=clusters[i]; 
      if (s<-1) printf("%d  %d %d \n",s,i,j);
      extract_row(cX,i,dp); scl_vec_mult(VE(Y,i),dp,dp); 
      vec_add(dp,cumA[j],cumA[j]); 

      if ((time==x[i])&(delta[i]==0))vec_add(qs,cumhatA[j],cumhatA[j]);  

      if (s<-1) print_vec(dp2); 
   }

   for (j=0;j<*antclust;j++) { 
      vec_add(cumhatA[j],cumA[j],dp1); 
      Mv(AI,dp1,dp2); replace_row(cumAt[j],s,dp2);  

      for(k=0;k<ps;k++) 
      for(c=0;c<ps;c++) ME(VAR,k,c)=ME(VAR,k,c)+VE(dp2,k)*VE(dp2,c); 

      if (*resample==1) 
      for (c=0;c<*px;c++) {l=j*(*px)+c; biid[l*(*Ntimes)+s]=VE(dp2,c);}
   }

   for (i=1;i<ps+1;i++) {
      var[i*(*Ntimes)+s]=ME(VAR,i-1,i-1); 
      est[i*(*Ntimes)+s]=VE(beta,i-1); score[i*(*Ntimes)+s]=VE(SCORE,i-1); }

} /* s=1 ... *Ntimes */ 


    if (*sim==1)
      comptestfunc(times,Ntimes,px,est,var,vcudif,antsim,test,testOBS,Ut,
		   simUt,cumAt,weighted,antclust,gamma2,line,timepowtest); 
  }
  // }}}
  else {
    // {{{ itfitsemi

    itfitsemi(times,
	      Ntimes,
	      x,
	      delta,
	      cause,
	      KMc,
	      z,
	      n,
	      px,
	      Nit,
	      score,
	      hess,
	      est,
	      var,
	      sim,
	      antsim,
	      test,
	      testOBS,
	      Ut,
	      simUt,
	      weighted,
	      gamma,
	      vargamma,
	      semi,
	      zsem,
	      pg,
	      trans,
	      gamma2,
	      CA,
	      line,
	      detail,
	      biid,
	      gamiid,
	      resample,
	      timepow,
	      clusters,
	      antclust,
	      timepowtest,
	      silent,
	      convc);
  }

  // }}}
    // {{{ convergence problems and free some stuff

  if (convproblems==1) silent[0]=2; 
  if (*semi==0) { 
    free_mats(&VAR,&X,&cX,&A,&AI,NULL); 
    if (*trans==2) {free_mats(&Z,NULL); free_vecs(&zgam,&gam,&zi,&rowZ,NULL);}

    free_vecs(&rr,&bhatub,&risk,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat,NULL); 
    free_vecs(&bet1,&xi,&rowX,NULL); 
    free_vecs(&dp,&dp1,&dp2,&dcovs,&pcovs,&betaub,&VdB,&qs,&SCORE,&beta,&difbeta,NULL); 

    for (i=0;i<*antclust;i++) {free_vec(cumhatA[i]); free_vec(cumA[i]); 
    free_mat(cumAt[i]);}
  }
free(vcudif);

// }}}
}
// {{{ itfitsemi
void itfitsemi(times,
	       Ntimes,
	       x,
	       delta,
	       cause,
	       KMc,
	       z,
	       antpers,
	       px,
	       Nit,
	       score,
	       hess,
	       est,
	       var,
	       sim,
	       antsim,
	       test,
	       testOBS,
	       Ut,
	       simUt,
	       weighted,
	       gamma,
	       vargamma,
	       semi,
	       zsem,
	       pg,
	       trans,
	       gamma2,
	       CA,
	       line,
	       detail,
	       biid,
	       gamiid,
	       resample,
	       timepow,
	       clusters,
	       antclust,
	       timepowtest,
	       silent,
	       convc)
double *times,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
*Ut,*simUt,*gamma,*zsem,*vargamma,*gamma2,*biid,*gamiid,*timepow,*timepowtest,
	*convc;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*weighted,
*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,*silent;
{
  // {{{ declare some stuff

  matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
  matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*dC,*XZ,*ZZ,*ZZI,*XZAI; 
  matrix *Ct,*C[*Ntimes],*Acorb[*Ntimes],*tmpM1,*tmpM2,*tmpM3,*tmpM4; 
  matrix *Vargam,*dVargam,*M1M2[*Ntimes],*Delta,*dM1M2,*M1M2t,*RobVargam;
  matrix *W3t[*antclust],*W4t[*antclust];
  vector *W2[*antclust],*W3[*antclust];
  vector *diag,*dB,*dN,*VdB,*AIXdN,*AIXlamt,*bhatt,*pbhat,*plamt;
  vector *korG,*pghat,*rowG,*gam,*dgam,*ZGdN,*IZGdN,*ZGlamt,*IZGlamt;
  vector *covsx,*covsz,*qs,*Y,*rr,*bhatub,*xi,*rowX,*rowZ,*difX,*zi,*z1,
    *tmpv1,*tmpv2,*lrisk;
  int sing,itt,i,j,k,l,s,c,pmax,totrisk,convproblems=0, 
      *n= calloc(1,sizeof(int)), *nx= calloc(1,sizeof(int)),
      *robust= calloc(1,sizeof(int));
  double time,dummy,dtime,timem;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double)),
	 *inc=calloc((*Ntimes)*(*px+1),sizeof(double));
  double lrr,fabs(), pow(); 
  int fixedcov; 
  // float gasdev(),expdev(),ran1();
  robust[0]=1; fixedcov=1; 
  n[0]=antpers[0]; nx[0]=antpers[0];
  timem=0; 

//if (*trans==1) for (j=0;j<*pg;j++) if (fabs(timepow[j]-1)>0.0001) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (fabs(timepow[j])>0.0001) {timem=1;break;}

  for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*px,W3t[j]);
    malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]);
    }

  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,&cdesignG,NULL); 
  malloc_mats(*px,*px,&tmpM1,&A,&AI,NULL);
  malloc_mats(*pg,*pg,&dVargam,&Vargam,&RobVargam,&tmpM2,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,&S,&ZZI,NULL); 
  malloc_mats(*px,*pg,&XZAI,&tmpM3,&Ct,&dC,&XZ,&dM1M2,&M1M2t,NULL);
  malloc_mat(*px,*pg,tmpM4); 
  for (j=0;j<*Ntimes;j++) { malloc_mat(*pg,*px,Acorb[j]); 
    malloc_mat(*px,*pg,C[j]); malloc_mat(*px,*pg,M1M2[j]);}
  malloc_mat(*Ntimes,*px,Delta); malloc_mat(*Ntimes,*px,tmpM1);

  malloc_vecs(*px,&covsx,&xi,&rowX,&difX,&tmpv1,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,NULL);
  malloc_vecs(*pg,&covsz,&zi,&rowZ,&tmpv2,&zi,&z1,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,&IZGlamt,NULL);
  malloc_vecs(*antpers,&Y,&bhatub,&rr,&lrisk,&dN,&pbhat,&pghat,&plamt,NULL);
  malloc_vec((*px)+(*pg),qs); 

  // }}}
			     
  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 

  if (fixedcov==1) {
    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*pg) ME(ldesignG,c,j)=zsem[j*(*antpers)+c]; } } }

  for (itt=0;itt<*Nit;itt++)
    {
      mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZGdN); vec_zeros(IZGlamt); 

      Mv(ldesignG,gam,pghat);
      for (s=0;s<*Ntimes;s++)
      {
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 

	  for(j=1;j<=*px;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];
	  Mv(ldesignX,bhatt,pbhat); 

	  totrisk=0; 
	  for (j=0;j<*antpers;j++) { 
	    VE(lrisk,j)=(x[j]>=time); totrisk=totrisk+VE(lrisk,j);
	    extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi); 

	    lrr=0; 
	    // {{{ compute P_1 and DP_1 
	    if (*trans==1) {
//	      if (timem>0)
		for (l=0;l<*pg;l++)
		  lrr=lrr+VE(gam,l)*VE(zi,l)*pow(time,timepow[l]); 
//	      else lrr=time*VE(pghat,j);   
	      VE(plamt,j)=1-exp(-VE(pbhat,j)-lrr); 
	      scl_vec_mult(1-VE(plamt,j),xi,xi);
	      scl_vec_mult((1-VE(plamt,j)),zi,zi);  
//	      if (timem>0) 
              for (l=0;l<*pg;l++) VE(zi,l)=pow(time,timepow[l])*VE(zi,l); 
//	      else scl_vec_mult(time,zi,zi); 
	    }
	    if (*trans==2) { 
//	      if (timem>0)  { 
		for (l=0;l<*pg;l++) lrr=lrr+VE(gam,l)*VE(zi,l)*pow(time,timepow[l]); 
//	    } 
//	      else lrr=VE(pghat,j);  
	      VE(rr,j)=exp(lrr);  
	      VE(plamt,j)=1-exp(-exp(VE(pbhat,j))*VE(rr,j)); 
	      scl_vec_mult((1-VE(plamt,j))*exp(VE(pbhat,j))*VE(rr,j),xi,xi); 
	      scl_vec_mult((1-VE(plamt,j))*exp(VE(pbhat,j))*VE(rr,j),zi,zi); 
//	      if (timem>0) 
	      for (l=0;l<*pg;l++) VE(zi,l)= pow(time,timepow[l])*VE(zi,l); 
	    }
            if (*trans==3) {
//              if (timem>0)  { 
	        for (l=0;l<*pg;l++) lrr=lrr+VE(gam,l)*VE(zi,l)*pow(time,timepow[l]); 
//	    } else lrr=VE(pghat,j);  
	      VE(rr,j)=exp(lrr);  
	      VE(plamt,j)=exp(VE(pbhat,j)+lrr)/(1+exp(VE(pbhat,j)+lrr)); 
	      dummy=VE(plamt,j)/(1+exp(VE(pbhat,j)+lrr)); 
	      scl_vec_mult(dummy,xi,xi); 
	      scl_vec_mult(dummy,zi,zi); 
//	      if (timem>0) 
   	    for (l=0;l<*pg;l++) VE(zi,l)= pow(time,timepow[l])*VE(zi,l); 
           }
	   if (*trans==4) {
//             if (timem>0)  { 
	        for (l=0;l<*pg;l++) lrr=lrr+VE(gam,l)*VE(zi,l)*pow(time,timepow[l]); 
//	   } else lrr=VE(pghat,j);  
	      VE(rr,j)=lrr;  
	      VE(plamt,j)=exp(VE(pbhat,j)+lrr); 
	      scl_vec_mult(VE(plamt,j),xi,xi); 
	      scl_vec_mult(VE(plamt,j),zi,zi); 
//	      if (timem>0) 
	      for (l=0;l<*pg;l++) VE(zi,l)= pow(time,timepow[l])*VE(zi,l); 
           }
	   // }}}

	    replace_row(cdesignX,j,xi); replace_row(cdesignG,j,zi); 
	    /*
	      if (itt==*Nit-1) {
	      if (KMc[j]<0.00001) vec_zeros(xi); else scl_vec_mult(1/KMc[j],xi,xi); 
	      scl_vec_mult(VE(lrisk,j),xi,xi); vec_add(xi,qs,qs); }
	    */
	    VE(Y,j)=((x[j]<=time) & (cause[j]==*CA))*1;
	    if (KMc[j]<0.001) VE(Y,j)=(VE(Y,j)/0.001)-VE(plamt,j); 
	    else VE(Y,j)=(VE(Y,j)/KMc[j])-VE(plamt,j);
	  }
	  MtA(cdesignX,cdesignX,A); 
	  invertS(A,AI,silent[0]); sing=0; 

          if (fabs(ME(AI,0,0))<.0000001) {
             convproblems=1; 
             if (*silent==0) printf("non-invertible design at time %lf\n",time); 
             itt=*Nit-1;  
	     for (k=1;k<=*px;k++) inc[k*(*Ntimes)+s]=0; 
          }

	  if (sing==0) { 
	  vM(cdesignX,Y,xi); Mv(AI,xi,AIXdN); 
	  MtA(cdesignG,cdesignG,ZZ); MtA(cdesignX,cdesignG,XZ);
	  MxA(AI,XZ,XZAI); MtA(XZAI,XZ,tmpM2); 
	  mat_subtr(ZZ,tmpM2,dCGam); 
	  scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

	  vM(cdesignG,Y,zi); vM(XZ,AIXdN,tmpv2); 
	  vec_subtr(zi,tmpv2,ZGdN); scl_vec_mult(dtime,ZGdN,ZGdN); 
	  vec_add(ZGdN,IZGdN,IZGdN); 
	  Acorb[s]=mat_transp(XZAI,Acorb[s]); 
	  C[s]=mat_copy(XZ,C[s]); 

	  /* scl_mat_mult(dtime,XZAI,tmpM4);mat_add(tmpM4,Ct,Ct); */
	  for (k=1;k<=*px;k++) inc[k*(*Ntimes)+s]=VE(AIXdN,k-1); 
	  }

	  if (itt==*Nit-1) {
	    for (i=0;i<*antpers;i++) 
            { // vec_zeros(tmpv1); vec_zeros(z1); 
              j=clusters[i]; 	
	      extract_row(cdesignX,i,xi); scl_vec_mult(VE(Y,i),xi,xi); 
	      Mv(AI,xi,rowX);
	      extract_row(cdesignG,i,zi); scl_vec_mult(VE(Y,i),zi,zi); 
	      vM(C[s],rowX,tmpv2); vec_subtr(zi,tmpv2,rowZ); 
	      scl_vec_mult(dtime,rowZ,rowZ); 
	     // vec_add(rowZ,z1,z1); 
	     // vec_add(rowX,tmpv1,tmpv1); 
	      vec_add(rowZ,W2[j],W2[j]); 
	      for (k=0;k<*px;k++) ME(W3t[j],s,k)= ME(W3t[j],s,k)+VE(rowX,k); 
	    }  
	 }
	} /* s=1,...Ntimes */

      invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZGdN,dgam); vec_add(gam,dgam,gam); 

      if (isnan(vec_sum(dgam)) && *silent==0) {
        if (convproblems==1) convproblems=3;  else convproblems=2; 
	printf("missing values in dgam %ld \n",(long int) s);
	vec_zeros(gam); 
      }

      dummy=0; for (k=0;k<*pg;k++)  dummy=dummy+fabs(VE(dgam,k)); 

      for (s=0;s<*Ntimes;s++) {
	vM(Acorb[s],dgam,korG); 
	est[s]=times[s]; var[s]=times[s]; 
	for (k=1;k<=*px;k++)  { est[k*(*Ntimes)+s]=
            est[k*(*Ntimes)+s]+inc[k*(*Ntimes)+s]-VE(korG,k-1); 
	  dummy=dummy+fabs(inc[k*(*Ntimes)+s]-VE(korG,k-1)); 
	  /* printf(" %lf ",est[k*(*Ntimes)+s]); printf(" \n");*/ }
      } /* s=1,...Ntimes */
      if (dummy<*convc && itt<*Nit-2) itt=*Nit-2; 

      if (*detail==1) { 
	printf(" iteration %d %d \n",itt,*Nit); 
	printf("Total score %lf \n",dummy); 
	printf(" gamma parmaeters \n"); print_vec(gam); 
	printf(" change in gamma \n"); print_vec(dgam); }

    } /*itt løkke */ 

  /* ROBUST VARIANCES   */ 
  if (*robust==1) 
    {
      for (s=0;s<*Ntimes;s++) {
	vec_zeros(VdB); 
	 for (i=0;i<*antclust;i++) {

	  Mv(ICGam,W2[i],tmpv2); vM(Acorb[s],tmpv2,rowX); 
	  extract_row(W3t[i],s,tmpv1); vec_subtr(tmpv1,rowX,difX); 
	  replace_row(W4t[i],s,difX); vec_star(difX,difX,tmpv1); 
	  vec_add(tmpv1,VdB,VdB);

	  if (*resample==1) {
	    if (s==1)
	      for (c=0;c<*pg;c++)
	      gamiid[c*(*antclust)+i]=gamiid[c*(*antclust)+i]+VE(tmpv2,c);
	    for (c=0;c<*px;c++) {l=i*(*px)+c; 
	      biid[l*(*Ntimes)+s]=biid[l*(*Ntimes)+s]+VE(difX,c);} }


	  if (s==0) { for (j=0;j<*pg;j++) for (k=0;k<*pg;k++) 
			ME(RobVargam,j,k)=ME(RobVargam,j,k)+VE(tmpv2,j)*VE(tmpv2,k);} 
	}  /* for (i=0;i<*antclust;i++) */ 
	for (k=1;k<*px+1;k++) var[k*(*Ntimes)+s]=VE(VdB,k-1); 

      } /* s=0..Ntimes*/
    }

  /* MxA(RobVargam,ICGam,tmpM2); MxA(ICGam,tmpM2,RobVargam);*/
  /* print_mat(RobVargam);  */ 

  for (j=0;j<*pg;j++) {gamma[j]=VE(gam,j);
    for (k=0;k<*pg;k++) {vargamma[k*(*pg)+j]=ME(RobVargam,j,k);}}
  
  if (convproblems==1) silent[0]=2; 
  if (*sim==1) {
    comptestfunc(times,Ntimes,px,est,var,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antclust,gamma2,line,timepowtest);
  }

  free_mats(&ldesignX,&A,&AI,&cdesignX,&ldesignG,&cdesignG,
	      &S,&dCGam,&CGam,&ICGam,&VarKorG,&dC,&XZ,&ZZ,&ZZI,&XZAI, 
	      &Ct,&tmpM1,&tmpM2,&tmpM3,&tmpM4,&Vargam,&dVargam,
	      &Delta,&dM1M2,&M1M2t,&RobVargam,NULL); 

  free_vecs(&qs,&Y,&rr,&bhatub,&diag,&dB,&dN,&VdB,&AIXdN,&AIXlamt,
	      &bhatt,&pbhat,&plamt,&korG,&pghat,&rowG,&gam,&dgam,&ZGdN,&IZGdN,
	      &ZGlamt,&IZGlamt,&xi,&rowX,&rowZ,&difX,&zi,&z1,&tmpv1,&tmpv2,&lrisk,
	      NULL); 

  for (j=0;j<*Ntimes;j++) {free_mat(Acorb[j]);free_mat(C[j]);free_mat(M1M2[j]);}
  for (j=0;j<*antclust;j++) {free_mat(W3t[j]); free_mat(W4t[j]);
    free_vec(W2[j]); free_vec(W3[j]); }
  free(vcudif); free(inc); 
  free(n); free(nx);  free(robust); 
}
// }}}

