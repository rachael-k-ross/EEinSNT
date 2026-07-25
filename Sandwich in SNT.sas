/********************************************************************/
/*------------------------------------------------------------------*/
/*																	*/
/* Illustration of sandwich variance estimation				        */
/* for sequential nested trials				                        */
/*																	*/
/* Rachael Ross July 2026											*/
/*																	*/
/*------------------------------------------------------------------*/
/********************************************************************/


/*##################################################-
#
# Data
#
##################################################-*/

/*/* dat_k1: baseline emulated trial (k=1) */*/
/*data dat_k1;*/
/*    call streaminit(7);*/
/*    do i = 1 to 100000;*/
/*        S  = 1;*/
/*        W  = rand('BERNOULLI', 0.2);*/
/*        A  = rand('BERNOULLI', 0.1 + 0.1*W);*/
/*        Y0 = rand('BERNOULLI', 0.3 + 0.1*W);*/
/*        Y1 = rand('BERNOULLI', 0.3 + 0.1*W - 0.1);*/
/*        Y  = A*Y1 + (1-A)*Y0;*/
/*        output;*/
/*    end;*/
/*run;*/
/**/
/*/* dat_k2: 2nd emulated trial, nested among A=0 from k1 */*/
/*data dat_k2;*/
/*    call streaminit(8);*/
/*    set dat_k1(where=(A=0) rename=(W=priorW));*/
/*    S  = 2;*/
/*    W  = rand('BERNOULLI', 0.1 + 0.4*priorW);*/
/*    A  = rand('BERNOULLI', 0.1 + 0.1*W);*/
/*    Y0 = rand('BERNOULLI', 0.3 + 0.1*W);*/
/*    Y1 = rand('BERNOULLI', 0.3 + 0.1*W - 0.1);*/
/*    Y  = A*Y1 + (1-A)*Y0;*/
/*    drop priorW;*/
/*run;*/
;

data dat_k1;
    input i W A Y;
    S = 1;
    datalines;
1  1 1 1
2  0 0 0
3  1 0 1
4  0 0 0
5  0 1 1
6  1 1 1
7  0 0 0
8  1 0 1
9  1 1 1
10 0 0 0
11 0 0 0
12 1 1 1
13 0 0 0
14 1 0 1
15 0 0 0
16 1 1 1
17 0 0 0
18 0 0 1
19 1 1 1
20 1 0 0
;
run;

data dat_k2;
    input i W A Y;
    S = 2;
    datalines;
2  0 0 1
3  1 1 1
4  0 0 0
7  1 0 1
8  0 1 1
10 0 0 0
11 1 0 0
13 0 0 0
14 1 1 1
15 0 0 1
17 1 1 1
18 0 0 0
20 1 0 0
;
run;

/* dat_long: stacked */
data dat_long;
    set dat_k1 dat_k2;
run;


/*##################################################-
#
# Sandwich with single trial
#
##################################################-*/


/*---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------*/

/* Save dataset to be used as df  */
data df;
set dat_k1; 
run;

/* Propensity score model  */
proc logistic data=df desc;              /* desc: model P(A=1) */
    model A = W;
    output out=ps_out predicted=ps;
    ods output ParameterEstimates=pe_ps;
run;

/* Calculate weights and risk difference */
data ipw;
    set ps_out;
    ipw = A/ps + (1-A)/(1-ps);
    contrib = ipw*Y*A - ipw*Y*(1-A);
run;

proc means data=ipw noprint;
    var contrib;
    output out=rd_out(keep=rd) mean=rd;
run;




/*---------------------------------------#
### Start IML for remaining steps	 ----
#---------------------------------------*/

proc iml;

	/*---------------------------------------#
	### Step 2 - Data setup for sandwich ----
	#---------------------------------------*/

    use df;
        read all var {A} into A;      /* trt vector */
        read all var {Y} into Y;      /* outcome vector */
        read all var {W} into W;
    close df;

    n = nrow(A);                      /* number of unique individuals */
    X = j(n,1,1) || W;                /* model matrix for ps model: intercept, W */


    /*---------------------------------------#
    ### Step 3 - Stack of ee's ----
    #---------------------------------------*/

    /* Estimating equation function: takes theta, returns n x p matrix */
    start eefx_single(theta) global(A, Y, X);
		theta = colvec(theta);						/* force theta to be column vector */

        alpha = theta[2:nrow(theta)];                          /* ps model parameters */

        pi    = 1 / (1 + exp(-(X*alpha)));            /* propensity score */
        ee_ps = (A - pi) # X;                         /* ee's for ps model */
        ee_rd = A#Y/pi - (1-A)#Y/(1-pi) - theta[1];   /* ee for risk difference */
		
        return (ee_rd || ee_ps);
    finish eefx_single;

    /* Column-sum wrapper, used for the bread (numerical derivative) */
    start sumcolee(theta) global(A, Y, X);
        ee = eefx_single(theta);
        return (ee[+,]);
    finish sumcolee;


    /*---------------------------------------#
    ### Step 4 - Sandwich estimator ----
    #---------------------------------------*/

    /* Read point estimates as vectors */
    use rd_out;
        read all var {rd} into rd_vec;      /* K x 1 */
    close rd_out;

    use pe_ps;
        read all var {Estimate} into alpha_vec;   /* (number of alpha params) x 1 */
    close pe_ps;

	theta_hat = rd_vec // alpha_vec; /* point estimates from Step 1 - 3 by 1 matrix */

    residuals = eefx_single(theta_hat);       /* residuals at the point estimates */
    meat      = (residuals` * residuals) / n; /* meat matrix */

    nparm = nrow(theta_hat);
    par   = nparm || . || .;
    call nlpfdd(func, deriv, na, "sumcolee", theta_hat, par);
    bread = -deriv / n;                       /* bread matrix */

    bread_inv = inv(bread);
    sandwich  = (bread_inv * meat * bread_inv`) / n;
    se        = sqrt(vecdiag(sandwich));

    /*---------------------------------------#
    ### Results ----
    #---------------------------------------*/

    b = theta_hat`;
    create out var {b se};
        append;
    close out;
quit;

data out;
    set out;
    lcl = b - 1.96*se;
    ucl = b + 1.96*se;
run;

proc print data=out; run;



/*##################################################-
#
# Sandwich with long dataset ----
# results still stratified by trial
#
##################################################-*/


/*---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------*/

/* Save dataset to be used as df  */
data df;	
set dat_long; 
run;

/* Propensity score model with W*factor(S) interaction, ref level S=1  */
proc logistic data=df desc;
    class S(param=ref ref='1');
    model A = W S W*S;
    output out=ps_out predicted=ps;
    ods output ParameterEstimates=pe_ps;
run;

proc print data=pe_ps; /* Confirm order of parameters - must match ee function below */
run;

/* Risk difference per trial */
data ipw;
    set ps_out;
    ipw = A/ps + (1-A)/(1-ps);
    contrib = ipw*Y*A - ipw*Y*(1-A);
run;

proc means data=ipw noprint nway;
    class S;
    var contrib;
    output out=rd_out(keep=S rd) mean=rd;
run;

proc sort data=rd_out; by S; run;   /* ensures rd_vec is ordered trial 1..K */


/* Sort data for sandwich implementation */
proc sort data=df; by i S; run;

/*---------------------------------------#
### Start IML for remaining steps	 ----
#---------------------------------------*/

proc iml;

	/*---------------------------------------#
	### Step 2 - Data setup for sandwich ----
	#---------------------------------------*/

   	use df;
        read all var {i S A W Y} into long;
    close df;

    i_long = long[,1];
    S_long = long[,2];
    A_long = long[,3];
    W_long = long[,4];
    Y_long = long[,5];

    n = ncol(unique(i_long));     /* number of unique individuals */
    K = ncol(unique(S_long));     /* number of trials */
	

	/* Model matrix for PS model: intercept, W, dummies for S=2..K, W*dummy for S=2..K */
    Xdum = j(nrow(long), K-1, 0);
    Xint = j(nrow(long), K-1, 0);
    do trial = 2 to K;
        Xdum[,trial-1] = (S_long = trial);
        Xint[,trial-1] = (S_long = trial) # W_long;
    end;
    X_long = j(nrow(long),1,1) || W_long || Xdum || Xint;


    /*---------------------------------------#
    ### Step 3 - Stack of ee's ----
    #---------------------------------------*/

	
    /* Estimating equation function: takes theta, returns n x p matrix */
	start eefx_snt(theta) global(A_long, Y_long, X_long, S_long, i_long, n, K);

		theta = colvec(theta);  /* Force to be column vector */

	    rd_theta = theta[1:K];
		alpha    = theta[(K+1):nrow(theta)];

	    pi_long = 1/(1+exp(-(X_long*alpha)));

	    long_ee_ps = (A_long - pi_long) # X_long;
	    long_ee_rd = A_long#Y_long/pi_long - (1-A_long)#Y_long/(1-pi_long);

	    ee_ps = j(n, ncol(X_long), 0);
	    ee_rd = j(n, K, 0);

	    do trial = 1 to K;                        
	        r = loc(S_long = trial);
	        p = i_long[r];
	        ee_ps[p, ] = ee_ps[p, ] + long_ee_ps[r, ];
	        ee_rd[p, trial] = long_ee_rd[r] - rd_theta[trial];
	    end;
		
	    return (ee_rd || ee_ps);
	finish eefx_snt;

    /* Column-sum wrapper, used for the bread (numerical derivative) */
    start sumcolee(theta) global(A_long, Y_long, X_long, S_long, i_long, n, K);
        ee = eefx_snt(theta);
        return (ee[+,]);
    finish sumcolee;

    /*---------------------------------------#
    ### Step 4 - Sandwich estimator ----
    #---------------------------------------*/

    /* Read point estimates as vectors */
    use rd_out;
        read all var {rd} into rd_vec;      /* K x 1 */
    close rd_out;

    use pe_ps;
        read all var {Estimate} into alpha_vec;   /* (number of alpha params) x 1 */
    close pe_ps;

	theta_hat = rd_vec // alpha_vec; /*  (K + p) x 1 column vector, matches theta ordering used in eefx_snt */
	
	/* Sandwich calcultaions */
    residuals = eefx_snt(theta_hat);       /* residuals at the point estimates */
    meat      = (residuals` * residuals) / n; /* meat matrix */

    nparm = nrow(theta_hat); print(nparm);
    par   = nparm || . || .;
    call nlpfdd(func, deriv, na, "sumcolee", theta_hat, par);
    bread = -deriv / n;                       /* bread matrix */

    bread_inv = inv(bread);
    sandwich  = (bread_inv * meat * bread_inv`) / n;
    se        = sqrt(vecdiag(sandwich));


    /*---------------------------------------#
    ### Results ----
    #---------------------------------------*/

    b = theta_hat`;
    create out var {b se};
        append;
    close out;
quit;

data out;
    set out;
    lcl = b - 1.96*se;
    ucl = b + 1.96*se;
run;

proc print data=out; run;



/*##################################################-
#
# Sandwich with pooling ----
#
##################################################-*/


/*---------------------------------------#
### Step 1 - Point estimates ----
#---------------------------------------*/

/* Save dataset to be used as df  */
data df;	
set dat_long; 
run;

/* Propensity score model, pooled -- W + factor(S), no interaction (matches model.matrix(~W + factor(S))) */
proc logistic data=df desc;
    class S(param=ref ref='1');
    model A = W S;
    output out=ps_out predicted=ps;
    ods output ParameterEstimates=pe_ps;
run;

proc print data=pe_ps; /* Confirm order of parameters - must match ee function below */
run;

/* Calculate weights and pooled risk difference (single rd across all person-trial rows) */
data ipw;
    set ps_out;
    ipw = A/ps + (1-A)/(1-ps);
    contrib = ipw*Y*A - ipw*Y*(1-A);
run;

proc means data=ipw noprint;
    var contrib;
    output out=rd_out(keep=rd) mean=rd;
run;

/* Sort data for sandwich implementation */
proc sort data=df; by i S; run;


/*---------------------------------------#
### Start IML for remaining steps	 ----
#---------------------------------------*/

proc iml;

	/*---------------------------------------#
	### Step 2 - Data setup for sandwich ----
	#---------------------------------------*/

   	use df;
        read all var {i S A W Y} into long;
    close df;

    i_long = long[,1];
    S_long = long[,2];
    A_long = long[,3];
    W_long = long[,4];
    Y_long = long[,5];

    n = ncol(unique(i_long));     /* number of unique individuals */
    K = ncol(unique(S_long));     /* number of trials */

    /* Model matrix for pooled PS model: intercept, W, dummies for S=2..K (no interaction) */
    Xdum = j(nrow(long), K-1, 0);
    do trial = 2 to K;
        Xdum[,trial-1] = (S_long = trial);
    end;
    X_long = j(nrow(long),1,1) || W_long || Xdum;


    /*---------------------------------------#
    ### Step 3 - Stack of ee's ----
    #---------------------------------------*/

	
 	/* Estimating equation function: takes theta, returns n x p matrix */
    start eefx_sntpool(theta) global(A_long, Y_long, X_long, S_long, i_long, n, K);

        theta = colvec(theta);              /* normalize orientation regardless of how nlpfdd passes it */

        rd     = theta[1];                  /* single pooled rd parameter */
        alpha  = theta[2:nrow(theta)];

        pi_long = 1/(1+exp(-(X_long*alpha)));

        long_ee_ps = (A_long - pi_long) # X_long;
        long_ee_rd = A_long#Y_long/pi_long - (1-A_long)#Y_long/(1-pi_long) - rd;

        ee_ps = j(n, ncol(X_long), 0);
        ee_rd = j(n, 1, 0);

        do trial = 1 to K;                  /* loop over trials, sum contributions within person */
            r = loc(S_long = trial);
            p = i_long[r];
            ee_ps[p, ]  = ee_ps[p, ] + long_ee_ps[r, ];
            ee_rd[p, 1] = ee_rd[p, 1] + long_ee_rd[r];   /* summed across trials -- pooled, not per-trial */
        end;

        return (ee_rd || ee_ps);
    finish eefx_sntpool;


    /* ALTERNATIVE eefx_sntpool that does not loop over K - this would be best for large K (e.g., >=50);
	NOTE: requires df sorted by i (ascending) */;
	/*start eefx_sntpool(theta) global(A_long, Y_long, X_long, S_long, i_long, n, K);

	    theta = colvec(theta);
	    rd    = theta[1];
	    alpha = theta[2:nrow(theta)];

	    pi_long = 1/(1+exp(-(X_long*alpha)));

	    long_ee_ps = (A_long - pi_long) # X_long;
	    long_ee_rd = A_long#Y_long/pi_long - (1-A_long)#Y_long/(1-pi_long) - rd;

	    n_long = nrow(i_long);

	    csum_rd = cusum(long_ee_rd);

	    p_ps = ncol(long_ee_ps);
	    csum_ps = j(n_long, p_ps, 0);
	    do col = 1 to p_ps;                     
	        csum_ps[,col] = cusum(long_ee_ps[,col]);
	    end;

	    is_end  = (i_long[2:n_long] ^= i_long[1:n_long-1]) // {1};
	    end_idx = colvec(loc(is_end));                     

	    n_end = nrow(end_idx);
	    prev_csum_rd = {0} // csum_rd[end_idx[1:n_end-1]];
	    prev_csum_ps = j(1, ncol(X_long), 0) // csum_ps[end_idx[1:n_end-1], ];

	    ee_rd = colvec(csum_rd[end_idx]) - colvec(prev_csum_rd);   
	    ee_ps = csum_ps[end_idx, ] - prev_csum_ps;

	    return (ee_rd || ee_ps);
	finish eefx_sntpool;*/


    /* Column-sum wrapper, used for the bread (numerical derivative) */
    start sumcolee(theta) global(A_long, Y_long, X_long, S_long, i_long, n, K);
        ee = eefx_sntpool(theta);
        return (ee[+,]);
    finish sumcolee;


    /*---------------------------------------#
    ### Step 4 - Sandwich estimator ----
    #---------------------------------------*/

 	use rd_out;
        read all var {rd} into rd_vec;      /* 1 x 1 */
    close rd_out;

    use pe_ps;
        read all var {Estimate} into alpha_vec;   /* (number of alpha params) x 1 */
    close pe_ps;

    theta_hat = rd_vec // alpha_vec; 
	

    residuals = eefx_sntpool(theta_hat);
    meat      = (residuals` * residuals) / n;

    nparm = nrow(theta_hat);
    par   = nparm || . || .;
    call nlpfdd(func, deriv, na, "sumcolee", theta_hat, par);
    bread = -deriv / n;

    bread_inv = inv(bread);
    sandwich  = (bread_inv * meat * bread_inv`) / n;
    se        = sqrt(vecdiag(sandwich));

    /*---------------------------------------#
    ### Results ----
    #---------------------------------------*/

    b = theta_hat`;
    create out var {b se};
        append;
    close out;
quit;

data out;
    set out;
    lcl = b - 1.96*se;
    ucl = b + 1.96*se;
run;

proc print data=out; run;
