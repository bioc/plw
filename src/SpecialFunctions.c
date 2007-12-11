#include<R.h>

// #include <stdio.h>
// #include <math.h>

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/

static const int my_incs[16] = {1073790977, 268460033, 67121153, 16783361, 4197377,
				1050113, 262913, 65921, 16577, 4193, 1073, 281, 77,
				23, 8, 1};


/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/

void getPosteriorMeansAndOdds(double *y,double *tau, double *p0, double *t1, double *t2, 
double *t3, double *m0, double *v0, int *p, int *n, int *k, int *k0,double *temp, double *post_mean, double *post_odds)
{
	int i,g;
	double kk,C,s0,temp0;
	
	kk=0.0-((*m0)+(*p))*0.5;

	for(g=0;g<(*n);g++){
		C=0.0;s0=0.0;
		for(i=0;i<(*k);i++){
			temp0=t1[g]-2.0*tau[i]*t2[g]+tau[i]*tau[i]*t3[0]+(*m0)*(*v0);
			temp[i]=pow(temp0,kk) *p0[i];
			C+=temp[i];
			s0+=(tau[i]*temp[i]);
		}
		post_mean[g]=s0/C;
		post_odds[g]=(C-temp[*k0])/temp[*k0];
	}
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void getPosteriorProbsSigmaHatAndMuhat(double *y,double *tau, double *p0, double *t1, double *t2, 
double *t3, double *m0, double *v0, int *p, int *n, int *k, int *k0,double *temp, double *n_hat, double *Sigma_hat, double *mu_hat)
{
	int i,j,g;
	double kk,C,s0,s1,s2,temp0,aa,sumw=0.0;
	double *yg;
	
	kk=0.0-((*m0)+(*p))*0.5;
	aa=(*m0)+(*p);
	mu_hat[0]=0.0;

	for(g=0;g<(*n);g++){
		C=0.0;s0=0.0;s1=0.0;s2=0.0;
		for(i=0;i<(*k);i++){
			temp0=t1[g]-2.0*tau[i]*t2[g]+tau[i]*tau[i]*t3[0]+(*m0)*(*v0);
			temp[i]=pow(temp0,kk) *p0[i];
			C+=temp[i];
			//temp0= E[I(gamma=tau_j)/c|y]*C
			temp0=aa*temp[i]/temp0;
			s0+=temp0;
			s1+=temp0*tau[i];
			s2+=temp0*tau[i]*tau[i];
		}
		s0=s0/C;
		s1=s1/C;
		s2=s2/C;
		//Vikten för att centrera : E[I(gamma=tau_k0)/c|y]
		temp0=t1[g]-2.0*tau[*k0]*t2[g]+tau[*k0]*tau[*k0]*t3[0]+(*m0)*(*v0);
		temp0=aa*temp[*k0]/(temp0*C);
		
		// Bidrag till n_hat:
		for(i=0;i<(*k);i++){
			n_hat[i]+=temp[i]/C;
		}
		yg=&y[g*(*p)];

		// Bidrag till mu_hat:
		mu_hat[0]+=(yg[(*p)-1]*temp0);
		sumw+=temp0;

		// Bidrag till Sigma_hat
		for(j=0;j<((*p)-1);j++){
			for(i=j;i<((*p)-1);i++){
				Sigma_hat[i*(*p)+j]+=(yg[i]*yg[j]*s0);
			}
			i=(*p)-1;
			Sigma_hat[i*(*p)+j]+=(yg[j]*(yg[i]*s0-s1));
		}
		i=(*p)-1;
		Sigma_hat[i*(*p)+i]+=(yg[i]*yg[i]*s0-2.0*yg[i]*s1+s2);
	}
	for(j=0;j<((*p)-1);j++){
		for(i=(j+1);i<(*p);i++){
			Sigma_hat[j*(*p)+i]=Sigma_hat[i*(*p)+j];
		}
	}
	mu_hat[0]=mu_hat[0]/sumw;
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void getPosteriorProbsAndSigmaHat(double *y,double *tau, double *p0, double *t1, double *t2, 
double *t3, double *m0, double *v0, int *p, int *n, int *k, int *k0,double *temp, double *n_hat, double *Sigma_hat)
{
	int i,j,g;
	double kk,C,s0,s1,s2,temp0,aa;
	double *yg;
	
	kk=0.0-((*m0)+(*p))*0.5;
	aa=(*m0)+(*p);

	for(g=0;g<(*n);g++){
		C=0.0;s0=0.0;s1=0.0;s2=0.0;
		for(i=0;i<(*k);i++){
			temp0=t1[g]-2.0*tau[i]*t2[g]+tau[i]*tau[i]*t3[0]+(*m0)*(*v0);
			temp[i]=pow(temp0,kk) *p0[i];
			C+=temp[i];
			//temp0= E[I(gamma=tau_j)/c|y]*C
			temp0=aa*temp[i]/temp0;
			s0+=temp0;
			s1+=temp0*tau[i];
			s2+=temp0*tau[i]*tau[i];
		}
		s0=s0/C;
		s1=s1/C;
		s2=s2/C;
		// Bidrag till n_hat:
		for(i=0;i<(*k);i++){
			n_hat[i]+=temp[i]/C;
		}
		// Bidrag till Sigma_hat
		yg=&y[g*(*p)];
		for(j=0;j<((*p)-1);j++){
			for(i=j;i<((*p)-1);i++){
				Sigma_hat[i*(*p)+j]+=(yg[i]*yg[j]*s0);
			}
			i=(*p)-1;
			Sigma_hat[i*(*p)+j]+=(yg[j]*(yg[i]*s0-s1));
		}
		i=(*p)-1;
		Sigma_hat[i*(*p)+i]+=(yg[i]*yg[i]*s0-2.0*yg[i]*s1+s2);
	}
	for(j=0;j<((*p)-1);j++){
		for(i=(j+1);i<(*p);i++){
			Sigma_hat[j*(*p)+i]=Sigma_hat[i*(*p)+j];
		}
	}
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   Version of above only using strata k0 to maximize Sigma
*/
void getPosteriorProbsAndSigmaHatStrata0(double *y,double *tau, double *p0, double *t1, double *t2, 
double *t3, double *m0, double *v0, int *p, int *n, int *k, int *k0,double *temp, double *n_hat, double *Sigma_hat)
{
	int i,j,g;
	double kk,C,s0,s1,s2,temp0,aa;
	double *yg;
	
	kk=0.0-((*m0)+(*p))*0.5;
	aa=(*m0)+(*p);

	for(g=0;g<(*n);g++){
		C=0.0;s0=0.0;s1=0.0;s2=0.0;
		for(i=0;i<(*k);i++){
			temp0=t1[g]-2.0*tau[i]*t2[g]+tau[i]*tau[i]*t3[0]+(*m0)*(*v0);
			temp[i]=pow(temp0,kk) *p0[i];
			C+=temp[i];
		}
		// Bidrag till n_hat:
		for(i=0;i<(*k);i++){
			n_hat[i]+=temp[i]/C;
		}
		// Bidrag till Sigma_hat
		temp0=aa*temp[*k0]/(C*(t1[g]+(*m0)*(*v0) ));
		yg=&y[g*(*p)];
		for(j=0;j<((*p));j++){
			for(i=j;i<((*p));i++){
				Sigma_hat[i*(*p)+j]+=(yg[i]*yg[j]*temp0);
			}
		}
	}
	for(j=0;j<((*p)-1);j++){
		for(i=(j+1);i<(*p);i++){
			Sigma_hat[j*(*p)+i]=Sigma_hat[i*(*p)+j];
		}
	}
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void getPosteriorProbs(double *tau, double *p0, double *temp1, double *temp2, double *temp3, double *m0, double *v0, int *p, int *n, int *k, int *k0,double *temp, double *p2, double *p3)
{
	int i,g;
	double kk,ss;
	kk=0.0-((*m0)+(*p))*0.5;

	for(g=0;g<(*n);g++){
		ss=0.0;
		for(i=0;i<(*k);i++){
			temp[i]=pow(temp1[g]-2.0*tau[i]*temp2[g]+tau[i]*tau[i]*temp3[0]+(*m0)*(*v0),kk) *p0[i];
			ss+=temp[i];
		}
		p3[g]=temp[(*k0)]/ss;
		for(i=0;i<(*k);i++){
			p2[i]+=temp[i]/ss;
		}
	}
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   covariance matrix mean equal zero.
   x a matrix (stored by col == by dimension) sized m*n
   sc a vector sized n
   covariance matrix sized m*m
*/
void cov_zero_mean_scaled_data(double *x,double *sc,int *n, int *m,double *covmat)
{
	int i,j,k;
	double sx=0.0,*xjn,*xkn;

	for(k=0;k<(*m);k++){
		for(j=0;j<=k;j++){
			sx=0.0;
			xjn=&x[j*(*n)];
			xkn=&x[k*(*n)];
			for(i=0;i<(*n);i++){
				sx=sx+xjn[i]*xkn[i]*sc[i];
			}
			covmat[k*(*m)+j]=sx/(double)(*n);
			covmat[j*(*m)+k]=sx/(double)(*n);
		}
	}	
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void getSS(double *A, double *y,int *n,int*p, double *res)
{
	int i,j,k;
	double ss=0.0,SS=0.0;
        double *yi; 

	/*Loop across rows of y */
	for(i=0;i<(*n);i++){
		yi=&y[i*(*p)];
		SS=0.0;
		for(j=0;j<(*p);j++){
			ss=0.0;
			for(k=0;k<(*p);k++)
				ss+=(A[j+k*(*p)]*yi[k]);
			SS+=(yi[j]*ss);
		}
		res[i]=SS;
	}
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void getSS_using_mu(double *A, double *mu, double *y,int *n,int*p, double *res)
{
	int i,j,k;
	double ss=0.0,SS=0.0;
        double *yi; 

	/*Loop across rows of y */
	for(i=0;i<(*n);i++){
		yi=&y[i*(*p)];
		SS=0.0;
		for(j=0;j<(*p);j++){
			ss=0.0;
			for(k=0;k<(*p);k++)
				ss+=(A[j+k*(*p)]*(yi[k]-mu[k]));
			SS+=((yi[j]-mu[j])*ss);
		}
		res[i]=SS;
	}
}
/* -------------------------------------------------------------------------------------------------------------------------------------------------
   Linear regression för varje rad i y.
   Antar att x har mean==0 och x^2 har summan 1
   x skall vara av längd n
   y av längd n*m
*/
void SimpLinReg(double *x, double *y,int *n, int *m,double *beta, double *t)
{
	int i,j;
	double sxy=0.0,sy=0.0,sy2=0.0;
        double *yi; 
	double s2;

	/*Loop across rows of yIterera*/
	for(i=0;i<(*m);i++){
		sxy=0.0;sy=0.0;sy2=0.0;
		yi=&y[i*(*n)];
		for(j=0;j<(*n);j++){
			sxy+=(x[j]*yi[j]);
			sy+=yi[j];
			sy2+=(yi[j]*yi[j]);
		}
		beta[i]=sxy;
		s2=(sy2-(sy*sy)/(double)(*n)-sxy*sxy)/(double)((*n)-2);
		t[i]=sxy*sxy/s2;
	}
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   mean and std for every row of x
   x a matrix (stored by row) sized n*m
   mean & sd vectors sized n
*/
void MeanAndSd(double *x,int *n, int *m,double *mean, double *sd2)
{
	int i,j;
	double sx=0.0,sx2=0.0;
	double *xi;

	/*Loop across rows of yIterera*/
	for(i=0;i<(*n);i++){
		sx=0.0;
		sx2=0.0;
		xi=&x[i*(*m)];
		for(j=0;j<(*m);j++){
			sx+=(xi[j]);
			sx2+=(xi[j]*xi[j]);
		}
		mean[i]=sx/(double)(*m);
		sd2[i]=(sx2-(sx*sx)/(double)(*m))/(double)((*m)-1);
	}
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   covariance matrix mean equal zero.
   x a matrix (stored by row == by observation) sized n*m
   covariance matrix sized m*m
*/
void CovMatrixZeroMean(double *x,int *n, int *m,double *covmat)
{
	int i,j,k;
	double sx=0.0;

	for(k=0;k<(*m);k++){
		for(j=0;j<=k;j++){
			sx=0.0;
			for(i=0;i<(*n);i++){
				sx=sx+x[i*(*m)+j]*x[i*(*m)+k];
			}
			covmat[k*(*m)+j]=sx/(double)(*n);
			covmat[j*(*m)+k]=sx/(double)(*n);
		}
	}	
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   Diag(A*B) where A and B are n*m and m*n respectively
   Matrixes stored by row
   Result is a vector sezed n
*/
void DiagAtimesB(double *a,double *b, int *n, int *m,double *diag)
{
	int i,j;
	double sx=0.0;

	for(i=0;i<(*n);i++){
		sx=0.0;
		for(j=0;j<(*m);j++){
			sx=sx+a[i*(*m)+j]*b[j*(*n)+i];
		}
		diag[i]=sx;
	}	
}


/* -------------------------------------------------------------------------------------------------------------------------------------------------
Gradient function for EstimateSigma.m.vbeta
*/
void grEstimateSigma(double *mat,double *delta, int *n, int *m,double *grad)
{
	int i,j;
	double sx=0.0, *mati;
	for(i=0;i<(*m);i++){
		sx=0.0;
		mati=&mat[i*(*n)];
		for(j=0;j<(*n);j++){
			sx=sx+mati[j]*delta[j];
		}
		grad[i]=sx/(double)(*n);
	}	
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   Diag(A*B) where A and B are n*m and m*n respectively
   Matrixes stored by row
   Result is a vector sezed n
*/
void DiagAtimesBv2(double *a,double *b, int *n, int *m,double *diag)
{
	int i,j;
	double sx=0.0,*ai,*bi;
	for(i=0;i<(*n);i++){
		sx=0.0;
		ai=&a[i*(*m)];
		bi=&b[i*(*m)];
		for(j=0;j<(*m);j++){
//			sx=sx+a[i*(*m)+j]*b[j+i*(*m)];
			sx=sx+ai[j]*bi[j];
		}
		diag[i]=sx;
	}	
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
   Linear regression för varje rad i y, OCH varje kolumns i x.
   Antar att varje kolumn i x har mean==0 och x^2 har summan 1
   x skall vara av längd n*p
   y av längd n*m
*/
void SimpLinReg2(double *x, double *y,int *n, int *m,int*p,double *beta, double *t,double *x_scale)
{
	int i,j,k;
	double sxy=0.0,sy=0.0,sy2=0.0;
        double *yi,*xi; 
	double s2;

	/*Loop across rows of yIterera*/
	for(k=0;k<(*p);k++){
		xi=&x[k*(*n)];
		for(i=0;i<(*m);i++){
			sxy=0.0;sy=0.0;sy2=0.0;
			yi=&y[i*(*n)];
			for(j=0;j<(*n);j++){
				sxy+=(xi[j]*yi[j]);
				sy+=yi[j];
				sy2+=(yi[j]*yi[j]);
			}
			beta[i+(*m)*k]=sxy/x_scale[k];
			s2=(sy2-(sy*sy)/(double)(*n)-sxy*sxy)/(double)((*n)-2);
			t[i+(*m)*k]=sxy*sxy/s2;
		}
	}
}

/*
int main()
{
	int n=11,m=2;
	double x[]={-0.5,-0.25,-0.25,-0.25,-0.25,0,0.25,0.25,0.25,0.25,0.5};
	double y[] ={0,1,2,3,4,5,6,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1};
	double beta[] = {0,0};
	double t[] = {0,0};

	SimpLinReg(x,y,&n,&m,beta, t);

	fprintf(stderr,"beta:  %f  %f\nt2:  %f  %f\n",beta[0],beta[1],t[0],t[1]);
	return 0;
}

*/




/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void my_sort(double *x, int *nn)
{
  double v;
  int i, j, h, t,n;

  n= (*nn);

  for (t = 0; my_incs[t] > n; t++);

  for (h = my_incs[t]; t < 16; h = my_incs[++t]){ 
    for (i = h; i < n; i++) { 
      v = x[i];
      j = i;
      while (j >= h && x[j - h] > v) {
	x[j] = x[j - h];
	j -= h;
      }
      x[j] = v;
    }
  }
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
double my_median(double *x,int *n)
{

  int index;

  if(n==0)
    return 0;

  index = (*n)/2;
  my_sort(x,n);

  if( index*2< (*n)){ 
    return x[index];
  }
  else{
    return ( ( x[index] + x[index-1] )*0.5 );
  }
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
double my_mad(double *x,int *n)
{
  double m;
  int i;
  
  m=my_median(x,n);
  
  for(i=0;i<(*n);i++)
     x[i]=fabs(x[i]-m);

  m=my_median(x,n);
  
  return(m);
}
	
/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
double my_mean(double *x,int *n)
{
  int i;
  double S=0.0;

  if(n==0)
    return 0;

  for(i=0;i<(*n);i++)
    S+=x[i];

  return ( S/(double)(*n));
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
double my_sd(double *x,int *n)
{
  int i;
  double S1=0.0,S2=0.0;

  if(n==0)
    return 0;

  for(i=0;i<(*n);i++){
    S1+=x[i];
	S2+=(x[i]*x[i]);
}
S1=sqrt((S2-S1*S1/(double)(*n))/(double)((*n)-1));

  return ( S1);
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void order_stat_by_index(double *x, int *index, int *q,int *n, int *length_res)
{
  int first=0,last=0,cur_index,cur_case=0,temp;
  double res;
  cur_index=index[first];

  while(first<(*n)){
    while((last<((*n)-1)) & (index[last+1]==cur_index))
      last++;
    temp=(last-first+1);
    res=my_median(&x[first],&temp);
    if(first+(*q)<(last+1)){
      res=x[first+(*q)];
    }
    else{
      res=x[last];
    }

    x[cur_case]=res;
    first=last+1;
    last=first;
    cur_index=index[first];
    cur_case++;
  }
  (*length_res)=cur_case;
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void median_by_index(double *x, int *index, int *n, int *length_res)
{
  int first=0,last=0,cur_index,cur_case=0,temp;
  double res;
  cur_index=index[first];

  while(first<(*n)){
    while((last<((*n)-1)) & (index[last+1]==cur_index))
      last++;
    temp=(last-first+1);
    res=my_median(&x[first],&temp);
    x[cur_case]=res;
    first=last+1;
    last=first;
    cur_index=index[first];
    cur_case++;
  }
  (*length_res)=cur_case;
}

/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void mad_by_index(double *x, int *index, int *n, int *length_res)
{
  int first=0,last=0,cur_index,cur_case=0,temp;
  double res;
  cur_index=index[first];

  while(first<(*n)){
    while((last<((*n)-1)) & (index[last+1]==cur_index))
      last++;
    temp=(last-first+1);
    res=my_mad(&x[first],&temp);
    x[cur_case]=res;
    first=last+1;
    last=first;
    cur_index=index[first];
    cur_case++;
  }
  (*length_res)=cur_case;
}


/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void mean_by_index(double *x, int *index, int *n, int *length_res)
{
  int first=0,last=0,cur_index,cur_case=0,temp;
  double res;
  cur_index=index[first];

  while(first<(*n)){
    while((last<((*n)-1)) & (index[last+1]==cur_index))
      last++;
    temp=(last-first+1);
    res=my_mean(&x[first],&temp);
    x[cur_case]=res;
    first=last+1;
    last=first;
    cur_index=index[first];
    cur_case++;
  }
  (*length_res)=cur_case;
}


/* -------------------------------------------------------------------------------------------------------------------------------------------------
*/
void sd_by_index(double *x, int *index, int *n, int *length_res)
{
  int first=0,last=0,cur_index,cur_case=0,temp;
  double res;
  cur_index=index[first];

  while(first<(*n)){
    while((last<((*n)-1)) & (index[last+1]==cur_index))
      last++;
    temp=(last-first+1);
    res=my_sd(&x[first],&temp);
    x[cur_case]=res;
    first=last+1;
    last=first;
    cur_index=index[first];
    cur_case++;
  }
  (*length_res)=cur_case;
}
