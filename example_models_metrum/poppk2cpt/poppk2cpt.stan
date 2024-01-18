functions{
  
  vector depot_2cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real q = params[2];
    real vc = params[3];
    real vp = params[4];
    real ka = params[5];
    
    real ke = cl/vc;
    real k_cp = q/vc;
    real k_pc = q/vp;
    
    vector[3] dydt;

    dydt[1] = -ka*y[1];                                // depot
    dydt[2] = ka*y[1] - (ke + k_cp)*y[2] + k_pc*y[3];  // central
    dydt[3] = k_cp*y[2] - k_pc*y[3];                   // peripheral
    
    return dydt;
  }
  
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> nSubjects;
  int nIIV;
  array[nObs] int<lower = 1> iObs;
  array[nSubjects] int<lower = 1> start;
  array[nSubjects] int<lower = 1> end;
  array[nt] int<lower = 1> cmt;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;
  vector<lower = 0>[nObs] cObs;
  
  int<lower = 1, upper = 3> solver;
  
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;
  int nCmt = 3;
  array[nSubjects] int nti;
  array[nCmt] real biovar;
  array[nCmt] real tlag;

  for(i in 1:nSubjects) nti[i] = end[i] - start[i] + 1;

  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0, upper = 500> CLHat;
  real<lower = 0, upper = 500> QHat;
  real<lower = 0, upper = 3500> V1Hat;
  real<lower = 0, upper = 3500> V2Hat;
  real<lower = 0, upper = 100> kaHat;
  real<lower = 0> sigma;
  
  // Inter-Individual variability
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0.01, upper = 2>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters{
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; // variable required for Matt's trick
  array[nTheta] real<lower = 0> theta;
  matrix<lower = 0>[nCmt, nt] x;
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;

  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = kaHat;

  // Matt's trick to use unit scale 
  thetaM = (rep_matrix(thetaHat, nSubjects) .* exp(diag_pre_multiply(omega, L * etaStd)))'; 
  
  for(j in 1:nSubjects)
  {
    theta[1] = thetaM[j, 1]; // CL
    theta[2] = thetaM[j, 2]; // Q
    theta[3] = thetaM[j, 3]; // V1
    theta[4] = thetaM[j, 4]; // V2
    theta[5] = thetaM[j, 5]; // ka

    if(solver == 1){
    
      x[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]], 
                                              amt[start[j]:end[j]],
                                              rate[start[j]:end[j]],
                                              ii[start[j]:end[j]],
                                              evid[start[j]:end[j]],
                                              cmt[start[j]:end[j]],
                                              addl[start[j]:end[j]],
                                              ss[start[j]:end[j]],
                                              theta, biovar, tlag);
    
    }else if(solver == 2){
      
      matrix[nCmt, nCmt] K = rep_matrix(0, nCmt, nCmt);
  
        K[1, 1] = -theta[5];
        K[2, 1] = theta[5];
        K[2, 2] = -((theta[1] + theta[2])/theta[3]);
        K[2, 3] = theta[2]/theta[4];
        K[3, 2] = theta[2]/theta[3];
        K[3, 3] = -theta[2]/theta[4];
        
        x[, start[j]:end[j]] =
          pmx_solve_linode(time[start[j]:end[j]],
                                amt[start[j]:end[j]],
                                rate[start[j]:end[j]],
                                ii[start[j]:end[j]],
                                evid[start[j]:end[j]],
                                cmt[start[j]:end[j]],
                                addl[start[j]:end[j]],
                                ss[start[j]:end[j]],
                                K, biovar, tlag);
      
    }else{
      
      x[, start[j]:end[j]] =
          pmx_solve_rk45(depot_2cmt_ode,
                         3,
                         time[start[j]:end[j]],
                         amt[start[j]:end[j]],
                         rate[start[j]:end[j]],
                         ii[start[j]:end[j]],
                         evid[start[j]:end[j]],
                         cmt[start[j]:end[j]],
                         addl[start[j]:end[j]],
                         ss[start[j]:end[j]],
                         theta, biovar, tlag);
      
    }
                                       
    cHat[start[j]:end[j]] = x[2, start[j]:end[j]] ./ theta[3]; // divide by V1
  }

  cHatObs  = cHat[iObs];
}

model{
  // Prior
  CLHat ~ lognormal(log(10), 0.25);
  QHat ~ lognormal(log(15), 0.5);
  V1Hat ~ lognormal(log(35), 0.25);
  V2Hat ~ lognormal(log(105), 0.5);
  kaHat ~ lognormal(log(2.5), 1);
  
  L ~ lkj_corr_cholesky(1);
  
  // Inter-individual variability (see transformed parameters block
  // for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);
  
  sigma ~ cauchy(0, 5);
  logCObs ~ normal(log(cHatObs), sigma);
}

