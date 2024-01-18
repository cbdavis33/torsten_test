// First Order Absorption (oral/subcutaneous)
// Two-compartment PK Model
// IIV on CL, VC, Q, VP, and Ka (full covariance matrix)
// proportional error - DV = IPRED*(1 + eps_p)
// Any of analytical, matrix-exponential, or general ODE solution (rk45 or bdf) 
//   using Torsten (user's choice)
// Implements threading for within-chain parallelization 
// Deals with BLOQ values by the "CDF trick" (M4)
// Since we have a normal distribution on the error, but the DV must be > 0, it
//   truncates the likelihood below at 0
// For PPC, it generates values from a normal that is truncated below at 0

functions{

  array[] int sequence(int start, int end) { 
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq; 
  } 
  
  int num_between(int lb, int ub, array[] int y){
    
    int n = 0;
    for(i in 1:num_elements(y)){
      if(y[i] >= lb && y[i] <= ub)
         n = n + 1;
    }
    return n;
    
  }
  
  array[] int find_between(int lb, int ub, array[] int y) {
    // vector[num_between(lb, ub, y)] result;
    array[num_between(lb, ub, y)] int result;
    int n = 1;
    for (i in 1:num_elements(y)) {
      if (y[i] >= lb && y[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  vector find_between_vec(int lb, int ub, array[] int idx, vector y) {
    
    vector[num_between(lb, ub, idx)] result;
    int n = 1;
    if(num_elements(idx) != num_elements(y)) reject("illegal input");
    for (i in 1:rows(y)) {
      if (idx[i] >= lb && idx[i] <= ub) {
        result[n] = y[i];
        n = n + 1;
      }
    }
    return result;
  }
  
  real normal_lb_rng(real mu, real sigma, real lb){
    
    real p_lb = normal_cdf(lb | mu, sigma);
    real u = uniform_rng(p_lb, 1);
    real y = mu + sigma * inv_Phi(u);
    return y; 

  }
  
  vector depot_2cmt_ode(real t, vector y, array[] real params, 
                        array[] real x_r, array[] int x_i){
    
    real cl = params[1];
    real vc = params[2];
    real q = params[3];
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
  
  real partial_sum_lpmf(array[] int seq_subj, int start, int end,
                        vector dv_obs, array[] int dv_obs_id, array[] int i_obs,
                        array[] real amt, array[] int cmt, array[] int evid, 
                        array[] real time, array[] real rate, array[] real ii, 
                        array[] int addl, array[] int ss,
                        array[] int subj_start, array[] int subj_end, 
                        vector CL, vector VC, vector Q, vector VP, vector KA, 
                        real sigma_p, 
                        vector lloq, array[] int bloq,
                        int n_random, int n_subjects, int n_total,
                        array[] real bioav, array[] real tlag, int n_cmt, 
                        int solver){
                           
    real ptarget = 0;
                              
    int N = end - start + 1;    // number of subjects in this slice  
    vector[n_total] dv_ipred;   
    matrix[n_total, n_cmt] x_ipred;
  
    int n_obs_slice = num_between(subj_start[start], subj_end[end], i_obs);
    array[n_obs_slice] int i_obs_slice = find_between(subj_start[start], 
                                                      subj_end[end], i_obs);
                                                
    vector[n_obs_slice] dv_obs_slice = find_between_vec(start, end, 
                                                        dv_obs_id, dv_obs);
    
    vector[n_obs_slice] ipred_slice;
    
    vector[n_obs_slice] lloq_slice = lloq[i_obs_slice];
    array[n_obs_slice] int bloq_slice = bloq[i_obs_slice];
    
    
    for(n in 1:N){            // loop over subjects in this slice
    
      int j = n + start - 1; // j is the ID of the current subject
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], Q[j], VC[j], VP[j], KA[j]})';
          
      }else if(solver == 2){
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
  
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -((CL[j] + Q[j])/VC[j]);
        K[2, 3] = Q[j]/VP[j];
        K[3, 2] = Q[j]/VC[j];
        K[3, 3] = -Q[j]/VP[j];
      
        x_ipred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
                           
      }else if(solver == 3){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_2cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q[j], VP[j], KA[j]}, bioav, tlag)';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(depot_2cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {CL[j], VC[j], Q[j], VP[j], KA[j]}, bioav, tlag)';
                         
      }
                      
      dv_ipred[subj_start[j]:subj_end[j]] = 
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
    
    }
  
    ipred_slice = dv_ipred[i_obs_slice];
    
    for(i in 1:n_obs_slice){
      real sigma_tmp = ipred_slice[i]*sigma_p;
      if(bloq_slice[i] == 1){
        ptarget += log_diff_exp(normal_lcdf(lloq_slice[i] | ipred_slice[i], 
                                                            sigma_tmp),
                                normal_lcdf(0.0 | ipred_slice[i], sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp); 
      }else{
        ptarget += normal_lpdf(dv_obs_slice[i] | ipred_slice[i], sigma_tmp) -
                   normal_lccdf(0.0 | ipred_slice[i], sigma_tmp);
      }
    }                                         
                              
    return ptarget;
                           
  }
  
}
data{
  
  int n_subjects;
  int n_total;
  int n_obs;
  array[n_obs] int i_obs;
  array[n_total] int ID;
  array[n_total] real amt;
  array[n_total] int cmt;
  array[n_total] int evid;
  array[n_total] real rate;
  array[n_total] real ii;
  array[n_total] int addl;
  array[n_total] int ss;
  array[n_total] real time;
  vector<lower = 0>[n_total] dv;
  array[n_subjects] int subj_start;
  array[n_subjects] int subj_end;
  vector[n_total] lloq;
  array[n_total] int bloq;
  
  real<lower = 0> location_tvcl;  // Prior Location parameter for CL
  real<lower = 0> location_tvvc;  // Prior Location parameter for VC
  real<lower = 0> location_tvq;   // Prior Location parameter for Q
  real<lower = 0> location_tvvp;  // Prior Location parameter for VP
  real<lower = 0> location_tvka;  // Prior Location parameter for KA
  
  real<lower = 0> scale_tvcl;     // Prior Scale parameter for CL
  real<lower = 0> scale_tvvc;     // Prior Scale parameter for VC
  real<lower = 0> scale_tvq;      // Prior Scale parameter for Q
  real<lower = 0> scale_tvvp;     // Prior Scale parameter for VP
  real<lower = 0> scale_tvka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_omega_cl; // Prior scale parameter for omega_cl
  real<lower = 0> scale_omega_vc; // Prior scale parameter for omega_vc
  real<lower = 0> scale_omega_q;  // Prior scale parameter for omega_q
  real<lower = 0> scale_omega_vp; // Prior scale parameter for omega_vp
  real<lower = 0> scale_omega_ka; // Prior scale parameter for omega_ka
  
  real<lower = 0> lkj_df_omega;   // Prior degrees of freedom for omega cor mat
  
  real<lower = 0> scale_sigma_p;  // Prior Scale parameter for proportional error
  
  int<lower = 0, upper = 1> prior_only; // Want to simulate from the prior?
  
  int<lower = 1, upper = 4> solver; // 1 = analytical, 2 = matrix exponential, 3 = rk45, 4 = bdf
 
}
transformed data{ 
  
  int grainsize = 1;
  
  vector<lower = 0>[n_obs] dv_obs = dv[i_obs];
  array[n_obs] int dv_obs_id = ID[i_obs];
  
  vector[n_obs] lloq_obs = lloq[i_obs];
  array[n_obs] int bloq_obs = bloq[i_obs];
  
  int n_random = 5;                    // Number of random effects
  int n_cmt = 3;                       // Number of states in the ODEs
  
  array[n_random] real scale_omega = {scale_omega_cl, scale_omega_vc, 
                                      scale_omega_q, scale_omega_vp,
                                      scale_omega_ka}; 
  
  array[n_subjects] int seq_subj = sequence(1, n_subjects); // reduce_sum over subjects
  
  array[n_cmt] real bioav = rep_array(1.0, n_cmt); // Hardcoding, but could be data or a parameter in another situation
  array[n_cmt] real tlag = rep_array(0.0, n_cmt);
  
}
parameters{ 
  
  real<lower = 0> TVCL;       
  real<lower = 0> TVVC; 
  real<lower = 0> TVQ;       
  real<lower = 0> TVVP;
  real<lower = 0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
    sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP))> TVKA;
  
  vector<lower = 0>[n_random] omega;
  cholesky_factor_corr[n_random] L;
  
  real<lower = 0> sigma_p;
  
  matrix[n_random, n_subjects] Z;
  
}
transformed parameters{
  
  vector[n_subjects] eta_cl;
  vector[n_subjects] eta_vc;
  vector[n_subjects] eta_q;
  vector[n_subjects] eta_vp;
  vector[n_subjects] eta_ka;
  vector[n_subjects] CL;
  vector[n_subjects] VC;
  vector[n_subjects] Q;
  vector[n_subjects] VP;
  vector[n_subjects] KA;

  {
  
    row_vector[n_random] typical_values = to_row_vector({TVCL, TVVC, TVQ, TVVP, 
                                                         TVKA});

    matrix[n_subjects, n_random] eta = diag_pre_multiply(omega, L * Z)';

    matrix[n_subjects, n_random] theta =
                          (rep_matrix(typical_values, n_subjects) .* exp(eta));
    
    eta_cl = col(eta, 1);
    eta_vc = col(eta, 2);
    eta_q = col(eta, 3);
    eta_vp = col(eta, 4);
    eta_ka = col(eta, 5);
    CL = col(theta, 1);
    VC = col(theta, 2);
    Q = col(theta, 3);
    VP = col(theta, 4);
    KA = col(theta, 5);
  
  }
  
}
model{ 
  
  // Priors
  TVCL ~ lognormal(log(location_tvcl), scale_tvcl);
  TVVC ~ lognormal(log(location_tvvc), scale_tvvc);
  TVQ ~ lognormal(log(location_tvq), scale_tvq);
  TVVP ~ lognormal(log(location_tvvp), scale_tvvp);
  TVKA ~ lognormal(log(location_tvka), scale_tvka) T[0.5*(TVCL/TVVC + TVQ/TVVC + TVQ/TVVP +
          sqrt((TVCL/TVVC + TVQ/TVVC + TVQ/TVVP)^2 - 4*TVCL/TVVC*TVQ/TVVP)), ];

  omega ~ normal(0, scale_omega);
  L ~ lkj_corr_cholesky(lkj_df_omega);
  
  sigma_p ~ normal(0, scale_sigma_p);
  
  to_vector(Z) ~ std_normal();
  
  // Likelihood
  if(prior_only == 0){
    target += reduce_sum(partial_sum_lupmf, seq_subj, grainsize,
                         dv_obs, dv_obs_id, i_obs,
                         amt, cmt, evid, time, 
                         rate, ii, addl, ss, subj_start, subj_end, 
                         CL, VC, Q, VP, KA,
                         sigma_p,
                         lloq, bloq,
                         n_random, n_subjects, n_total,
                         bioav, tlag, n_cmt, solver);
  }
}
generated quantities{
  
  real<lower = 0> sigma_sq_p = square(sigma_p);

  real<lower = 0> omega_cl = omega[1];
  real<lower = 0> omega_vc = omega[2];
  real<lower = 0> omega_q = omega[3];
  real<lower = 0> omega_vp = omega[4];
  real<lower = 0> omega_ka = omega[5];

  real<lower = 0> omega_sq_cl = square(omega_cl);
  real<lower = 0> omega_sq_vc = square(omega_vc);
  real<lower = 0> omega_sq_q = square(omega_q);
  real<lower = 0> omega_sq_vp = square(omega_vp);
  real<lower = 0> omega_sq_ka = square(omega_ka);

  real cor_cl_vc;
  real cor_cl_q;
  real cor_cl_vp;
  real cor_cl_ka;
  real cor_vc_q;
  real cor_vc_vp;
  real cor_vc_ka;
  real cor_q_vp;
  real cor_q_ka;
  real cor_vp_ka;
  real omega_cl_vc;
  real omega_cl_q;
  real omega_cl_vp;
  real omega_cl_ka;
  real omega_vc_q;
  real omega_vc_vp;
  real omega_vc_ka;
  real omega_q_vp;
  real omega_q_ka;
  real omega_vp_ka;

  vector[n_obs] ipred;
  vector[n_obs] pred;
  vector[n_obs] dv_ppc;
  vector[n_obs] log_lik;
  vector[n_obs] res;
  vector[n_obs] wres;
  vector[n_obs] ires;
  vector[n_obs] iwres;
 
  {

    matrix[n_random, n_random] R = multiply_lower_tri_self_transpose(L);
    matrix[n_random, n_random] Omega = quad_form_diag(R, omega);

    vector[n_total] dv_pred;
    matrix[n_total, n_cmt] x_pred;
    vector[n_total] dv_ipred;
    matrix[n_total, n_cmt] x_ipred;

    cor_cl_vc = R[1, 2];
    cor_cl_q = R[1, 3];
    cor_cl_vp = R[1, 4];
    cor_cl_ka = R[1, 5];
    cor_vc_q = R[2, 3];
    cor_vc_vp = R[2, 4];
    cor_vc_ka = R[2, 5];
    cor_q_vp = R[3, 4];
    cor_q_ka = R[3, 5];
    cor_vp_ka = R[4, 5];

    omega_cl_vc = Omega[1, 2];
    omega_cl_q = Omega[1, 3];
    omega_cl_vp = Omega[1, 4];
    omega_cl_ka = Omega[1, 5];
    omega_vc_q = Omega[2, 3];
    omega_vc_vp = Omega[2, 4];
    omega_vc_ka = Omega[2, 5];
    omega_q_vp = Omega[3, 4];
    omega_q_ka = Omega[3, 5];
    omega_vp_ka = Omega[4, 5];

    for(j in 1:n_subjects){
      
      if(solver == 1){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {CL[j], Q[j], VC[j], VP[j], KA[j]})';
                           
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_twocpt(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           {TVCL, TVQ, TVVC, TVVP, TVKA})';
          
      }else if(solver == 2){
        
        real tvke = TVCL/TVVC;
        real tvk_cp = TVQ/TVVC;
        real tvk_pc = TVQ/TVVP;
        
        matrix[n_cmt, n_cmt] K = rep_matrix(0, n_cmt, n_cmt);
        matrix[n_cmt, n_cmt] K_tv = rep_matrix(0, n_cmt, n_cmt);
        
        K[1, 1] = -KA[j];
        K[2, 1] = KA[j];
        K[2, 2] = -((CL[j] + Q[j])/VC[j]);
        K[2, 3] = Q[j]/VP[j];
        K[3, 2] = Q[j]/VC[j];
        K[3, 3] = -Q[j]/VP[j];
      
        x_ipred[subj_start[j]:subj_end[j], ] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K, bioav, tlag)';
        
        K_tv[1, 1] = -TVKA;
        K_tv[2, 1] = TVKA;
        K_tv[2, 2] = -((TVCL + TVQ)/TVVC);
        K_tv[2, 3] = TVQ/TVVP;
        K_tv[3, 2] = TVQ/TVVC;
        K_tv[3, 3] = -TVQ/TVVP;

        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_linode(time[subj_start[j]:subj_end[j]],
                           amt[subj_start[j]:subj_end[j]],
                           rate[subj_start[j]:subj_end[j]],
                           ii[subj_start[j]:subj_end[j]],
                           evid[subj_start[j]:subj_end[j]],
                           cmt[subj_start[j]:subj_end[j]],
                           addl[subj_start[j]:subj_end[j]],
                           ss[subj_start[j]:subj_end[j]],
                           K_tv, bioav, tlag)';
                           
      }else if(solver == 3){
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_2cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {CL[j], VC[j], Q[j], VP[j], KA[j]}, 
                         bioav, tlag)';
                         
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_rk45(depot_2cmt_ode,
                         n_cmt,
                         time[subj_start[j]:subj_end[j]],
                         amt[subj_start[j]:subj_end[j]],
                         rate[subj_start[j]:subj_end[j]],
                         ii[subj_start[j]:subj_end[j]],
                         evid[subj_start[j]:subj_end[j]],
                         cmt[subj_start[j]:subj_end[j]],
                         addl[subj_start[j]:subj_end[j]],
                         ss[subj_start[j]:subj_end[j]],
                         {TVCL, TVVC, TVQ, TVVP, TVKA}, 
                         bioav, tlag)';
                           
      }else{
        
        x_ipred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(depot_2cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {CL[j], VC[j], Q[j], VP[j], KA[j]}, 
                        bioav, tlag)';
                        
        x_pred[subj_start[j]:subj_end[j],] =
          pmx_solve_bdf(depot_2cmt_ode,
                        n_cmt,
                        time[subj_start[j]:subj_end[j]],
                        amt[subj_start[j]:subj_end[j]],
                        rate[subj_start[j]:subj_end[j]],
                        ii[subj_start[j]:subj_end[j]],
                        evid[subj_start[j]:subj_end[j]],
                        cmt[subj_start[j]:subj_end[j]],
                        addl[subj_start[j]:subj_end[j]],
                        ss[subj_start[j]:subj_end[j]],
                        {TVCL, TVVC, TVQ, TVVP, TVKA}, 
                        bioav, tlag)';                
                         
      }
      
      dv_ipred[subj_start[j]:subj_end[j]] =
        x_ipred[subj_start[j]:subj_end[j], 2] ./ VC[j];
      
      dv_pred[subj_start[j]:subj_end[j]] =
        x_pred[subj_start[j]:subj_end[j], 2] ./ TVVC;
      
    }

    pred = dv_pred[i_obs];
    ipred = dv_ipred[i_obs];

  }

  res = dv_obs - pred;
  ires = dv_obs - ipred;

  for(i in 1:n_obs){
    real ipred_tmp = ipred[i];
    real sigma_tmp = ipred_tmp*sigma_p;
    dv_ppc[i] = normal_lb_rng(ipred_tmp, sigma_tmp, 0.0);
    if(bloq_obs[i] == 1){
      // log_lik[i] = log(normal_cdf(lloq_obs[i] | ipred_tmp, sigma_tmp) -
      //                  normal_cdf(0.0 | ipred_tmp, sigma_tmp)) -
      //              normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
      log_lik[i] = log_diff_exp(normal_lcdf(lloq_obs[i] | ipred_tmp, sigma_tmp),
                                normal_lcdf(0.0 | ipred_tmp, sigma_tmp)) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
    }else{
      log_lik[i] = normal_lpdf(dv_obs[i] | ipred_tmp, sigma_tmp) -
                   normal_lccdf(0.0 | ipred_tmp, sigma_tmp);
    }
    wres[i] = res[i]/sigma_tmp;
    iwres[i] = ires[i]/sigma_tmp;
  }
  
}

