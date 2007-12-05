#include<R.h>


static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}

static double ** matrix(int nrh, int nch)
{
    int   i;
    double **m;

    m = (double **) R_alloc((nrh + 1), sizeof(double *));
    for (i = 0; i <= nrh; i++)
	m[i] = (double*) R_alloc((nch + 1), sizeof(double));
    return m;
}

static double ** Lmatrix(int n)
{
    int   i;
    double **m;

    m = (double **) R_alloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
	m[i] = (double *) R_alloc((i + 1), sizeof(double));
    return m;
}

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0

void  functionANDgradient(int n0, double *b, double *mat2, double *delta, double *meanmat2,int nn, double *f, double *gr,double *bb)
{
	int i,j;
	double temp,temp2;
	*f=0.0;
	for(j=0;j<n0;j++){	gr[j]=0.0; bb[j]=b[j];	}

	for(i=0;i<nn;i++){
		temp=0.0;
		for(j=0;j<n0;j++)		temp+=(mat2[i+j*nn]*b[j]);
		temp2=exp(temp)*delta[i];
		f[0]+=(temp-temp2);
		for(j=0;j<n0;j++)		gr[j]+=(mat2[i+j*nn]*temp2);		
	}
	f[0]=0-f[0]/(double)(nn);
	for(j=0;j<n0;j++)		gr[j]=	gr[j]/(double)(nn) - meanmat2[j];

}

int  vectorsequal(int n0, double *a, double *b)
{
	int i,res=1;
	for(i=0;i<n0;i++)
		if(a[i]!=b[i])
			res=res*0;
	return (res);
}
		
	

/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */



void
SpecialOptim(int *n0, double *b, double *Fmin, double *mat2, double *delta, double *meanmat2,int *nn,
      int *maxit, int *trace,
      double *abstol, double *reltol,
      int *fncount, int *grcount, int *fail)
{
    int accpoint, enough;
    double *g, *t, *X, *c, **B,*bf,*bg;
    int   count, funcount, gradcount;
    double f, gradproj,fg;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n;

    if (maxit <= 0) {
	*fail = 0;
//	*Fmin = fminfn(n0, b, ex);
	*fncount = *grcount = 0;
	return;
    }


    n = *n0;
    g = (double *)vect(n);
    t = vect(n);
    bf = vect(n);
    bg = vect(n);
    X = vect(n);
    c = vect(n);
    B = Lmatrix(n);
//    f = fminfn(n0, b, ex);   //******************** f
    functionANDgradient(n,b,mat2,delta,meanmat2,*nn,&f,g,bf);
    for(i=0;i<n;i++)	bg[i]=bf[i];
//    if (!R_FINITE(f))
//	fprintf(stderr,"initial value in 'vmmin' is not finite\n");
    if (*trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
//    fmingr(n0, b, g, ex);   //******************** GR
// f=*Fmin=f(b), g=gr(b)
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[i];
	    c[i] = g[i];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[j];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[j];
	    t[i] = s;
	    gradproj += s * g[i];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = 0;
	    do {
		count = 0;
		for (i = 0; i < n; i++) {
		    b[i] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[i]) /* no change */
			count++;
		}
		if (count < n) {
//		    f = fminfn(n0, b, ex);  //******************** f
		    functionANDgradient(n,b,mat2,delta,meanmat2,*nn,&f,g,bf);
		    funcount++;
		    accpoint = //R_FINITE(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (accpoint==0) {
			steplength *= stepredn;
		    }
		}
	    } while (!(count == n || accpoint==1));
	    enough = (f > (*abstol)) &&
		fabs(f - *Fmin) > (*reltol) * (fabs(*Fmin) + (*reltol));
	    /* stop if value if small or if relative change is low */
	    if (enough==0) {
		count = n;
		*Fmin = f;
	    }
	    if (count < n) {/* making progress */
		*Fmin = f;
//		fmingr(n0, b, g, ex);     //******************** GR
		if(vectorsequal(n,b,bg)==0){
		    functionANDgradient(n,b,mat2,delta,meanmat2,*nn,&fg,g,bg);

	   	    gradcount++;
		}
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[i] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (*trace && (iter % 10 == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= (*maxit)) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
    if (*trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < (*maxit)) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < (*maxit)) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
}
