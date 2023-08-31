// 
// statistical, seperable age-structured population model for GOA dusky rockfish
// based upon previous "base" model though has 1,150 fewer lines
// ben.williams@noaa.gov
// 2023-08

DATA_SECTION

  !!CLASS ofstream evalout("evalout.prj");

    // read the control file
    !! ad_comm::change_datafile_name("goa_dr_2022.ctl");    
    !! *(ad_comm::global_datafile) >>  model_name; 
    !! *(ad_comm::global_datafile) >>  data_file; 

    // phases - general
    init_int    ph_Fdev               // fishing mortality deviations
    init_int    ph_avg_F              // estimating average fishing mortality
    init_int    ph_recdev             // estimating recruitment deviations
    init_int    ph_fish_sel           // estimating fishing selectivity
    init_int    ph_srv1_sel           // estimating survey selectivity
    init_int    ph_F50                // last phase run F50 (4)

    // priors, cvs, phases for specific parameters
    init_number mprior                // prior mean for natural mortality
    init_number cvmprior              // prior CV for natural mortality
    init_int    ph_m                  // phase for estimating natural mortality

    init_int    ph_Rzero              // phase for estimating Rzero (determined by model specification

    init_number sigrprior             // prior mean for recruitment deviations
    init_number cvsigrprior           // prior CV for recruitment deviations
    init_int    ph_sigr               // phase for recruiment deviations

    init_number q_srv1prior           // prior mean for catchability coefficient
    init_number cvq_srv1prior         // prior CV for catchability coefficient
    init_int    ph_q_srv1             // phase for estimating catchability
    
    // weights
    init_int    yr_catchwt            // year catch-wt changes - early years are suspect.
    init_number wt_ssqcatch           // catch estimation - early 
    init_number wt_ssqcatch2          // catch estimation - later
    init_number wt_srv1               // survey biomass estimation
    
    init_number wt_fish_age           // fishery age compositions
    init_number wt_srv1_age           // survey age compositions
    init_number wt_fish_size          // fishery size compositions
    init_number wt_srv1_size          // survey size compositions

    init_number wt_rec_var            // Weight for estimation recruitment variations penalty
    init_number wt_sel_reg_fish       // Weight on fishery selectivity regularity penalty
    init_number wt_sel_reg_srv1       // Weight on survey selectivity regularity penalty

    init_number wt_fmort_reg          // Weight on fishing mortality regularity
    init_number wt_avg_sel            // Average selectivity penalty
    init_number initial_LMR           // weight of penalty to make mean recruitment like Rzerox
    init_number wt_Rzero              // weight of penalty to make mean recruitment like Rzerox
    init_number yieldratio  

    // transform
    number log_mprior
    !! log_mprior = log(mprior);
    number log_q_srv1prior       
    !! log_q_srv1prior=log(q_srv1prior);

    // ----------------------------------------------
    // read from .dat file
        !! ad_comm::change_datafile_name(data_file);    // Read data from the data file

        // Start and end years, recruitment age, number of age and length bins
        init_int      styr
        init_int      endyr
        init_int      recage
        init_int      nages_D
        init_int      nages_M
        init_int      nlenbins
        init_int      n_ageage_mat
        init_int      n_sizeage_mat
        init_vector   len_bin_labels(1,nlenbins)
        init_number spawn_fract; // Spawning Month
        init_vector   p_mature(1,nages_M) // maturity at age
        init_vector   wt(1,nages_M) // weight at age 

        // Observed catches
        init_vector   obs_catch_early(styr,yr_catchwt)
        init_vector   obs_catch_later(yr_catchwt+1,endyr)

        // Trawl Survey biomass estimates
        init_int      nyrs_srv1                     // number of years of survey biomass estimates
        init_ivector  yrs_srv1(1,nyrs_srv1)         // years survey conducted in
        init_vector   obs_srv1_biom(1,nyrs_srv1)    // mean estimate of biomass
        init_vector   obs_srv1_se(1,nyrs_srv1)      // standard error of survey biomass estimates
        init_vector   obs_srv1_lci(1,nyrs_srv1)     // lower confidence interval, for graphing not used in estimation
        init_vector   obs_srv1_uci(1,nyrs_srv1)     // upper confidence interval
        
        // Fishery age composition data
        init_int             nyrs_fish_age                           // number of years of fishery age compos
        init_ivector         yrs_fish_age(1,nyrs_fish_age)           // the years of age comps
        init_vector          nsamples_fish_age(1,nyrs_fish_age)      // total sample size for each age comp.
        init_vector          nhauls_fish_age(1,nyrs_fish_age)        // number of hauls for for each age comp.
        init_vector          nmulti_fish_age(1,nyrs_fish_age)        // input ESS for fishery age
        init_ivector         age_age_ind_fsh(1,nyrs_fish_age)        // some measure of relative sample size for each fishery comp
        init_matrix          oac_fish(1,nyrs_fish_age,1,nages_D)     // the actual year by year age comps
        //  vector               nmulti_fish_age(1,nyrs_fish_age)        // Input N for fishery age

        // Bottom trawl survey age composition data
        init_int             nyrs_srv1_age                           // number of years of survey age compositions
        init_ivector         yrs_srv1_age(1,nyrs_srv1_age)           // the years of survey age comps
        init_vector          nsamples_srv1_age(1,nyrs_srv1_age)      // total sample size for each age comp.
        init_vector          nhauls_srv1_age(1,nyrs_srv1_age)        // number of hauls for for each age comp.
        init_vector          nmulti_srv1_age(1,nyrs_srv1_age)        // input ESS for survey age
        init_ivector         age_age_ind_srv(1,nyrs_srv1_age)        // some measure of relative sample size for each fishery comp
        init_matrix          oac_srv1(1,nyrs_srv1_age,1,nages_D)     // the year by year age survey age comps
        //  vector               nmulti_srv1_age(1,nyrs_srv1_age)        // Input N for survey 1 age

        // Fishery size composition data
        init_int             nyrs_fish_size                          // number of years of fishery size comps
        init_ivector         yrs_fish_size(1,nyrs_fish_size)         // the years of fishery size comps
        init_vector          nsamples_fish_size(1,nyrs_fish_size)    // totals for fish lengths
        init_vector          nhauls_fish_size(1,nyrs_fish_size)      // hauls for fish lengths by year
        init_vector          nmulti_fish_size(1,nyrs_fish_size)      // input ESS for fishery size
        init_ivector         siz_age_ind_fsh(1,nyrs_fish_size)       // some measure of relative sample size for each fishery comp
        init_matrix          osc_fish(1,nyrs_fish_size,1,nlenbins)   // year by year fishery size comps
        //  vector               nmulti_fish_size(1,nyrs_fish_size)      // Input N for fishery size

        // Bottom trawl survey size composition data
        init_int             nyrs_srv1_size                          // number of years of survey size comps
        init_ivector         yrs_srv1_size(1,nyrs_srv1_size)         // the years of survey size comps
        init_vector          nsamples_srv1_size(1,nyrs_srv1_size)    // total lengths for survey 1 by year
        init_vector          nhauls_srv1_size(1,nyrs_srv1_size)      // total hauls for length samples in survey 1 by year
        init_ivector         siz_age_ind_srv1(1,nyrs_srv1_size)      // some measure of relative sample size for each fishery comp
        init_matrix          osc_srv1(1,nyrs_srv1_size,1,nlenbins)   // year by year size comps

        // Size-age transition matrix:  proportion at size given age
        init_3darray  sizeage(1,n_sizeage_mat,1,nages_M,1,nlenbins)       //size comp #1

        // Ageing error transition matrix:  proportion at reader age given true age
        init_3darray  ageage(1,n_ageage_mat,1,nages_M,1,nages_D)            // ageing error matrix  

        //end of .dat file marker
        init_int eof;

        // transform & containers
        int styr_rec
        int styr_sp
        int nyrs
        !!  nyrs = endyr - styr + 1;
        !! styr_rec = (styr - nages_M) + 1;     // First year of recruitment
        !! styr_sp  = styr_rec - recage ;     // First year of spawning biomass  

        !! spawn_fract = (spawn_fract - 1) / 12;
        vector wt_mature(1,nages_M);                  // Weight of mature fish vector at age
        
        // Declare index variables
        int                  iyr
        int                  i
        int                  j
        int                  l;
        vector yy(styr,endyr);
        !! yy.fill_seqadd(styr,1) ;
        vector aa(1,nages_M);
        !! aa.fill_seqadd(recage,1) ;
        vector aa_D(1,nages_D);
        !! aa_D.fill_seqadd(recage,1) ;

    // ========================================================================
    // checking that data file is correct
 LOCAL_CALCS

   if(eof==42) cout << "The data has been read correctly!";
   else { cout << "You f'ed up your data file!" << "\n"; exit(1); }
   if(wt_rec_var==0) {
     if (ph_sigr>0) {
       cout << "Warning, wt_rec_var is zero, so can't estimate sigr!@" << "\n";
       cout << "turning sigr off " << "\n";
       ph_sigr =-4;
       cout << "hit any key, then enter to continue" << "\n";
       char  xxx; 
       cin >> xxx;
     }
   }

  for (int i=1;i<=nages_M;i++){
  wt_mature(i) = wt(i) * p_mature(i)/2;
  }
 END_CALCS

INITIALIZATION_SECTION

    // Starting values for estimated parameters; these values over-ride all others
        logm                    log_mprior
        log_mean_rec            initial_LMR
        sigr                    sigrprior    
        a50                     7.5
        delta                   3.
        a50_srv1                7.3
        delta_srv1              3.8

PARAMETER_SECTION
    // Stock-recruitment
    vector               Sp_Biom(styr,endyr)
    init_number          log_mean_rec(1);                        // Unfish equil recruitment (logged)
    init_bounded_number  sigr(0.3,10,ph_sigr);                   // Recruitment sderr parameter

    // Fishery selectivity
    init_number          a50(ph_fish_sel);                       // age at 50% selection                                                   
    init_number          delta(ph_fish_sel);                     // age between 50% selection and 95% selection....
    vector               fish_sel(1,nages_M);                    // vectory of fishery selectivty parameters on arithmetic scale

    // Trawl Survey selectivity
    init_number          a50_srv1(ph_srv1_sel);                  // age at 50% selection                                                   
    init_number          delta_srv1(ph_srv1_sel);                // age between 50% selection and 95% selection....
    vector               srv1_sel(1,nages_M);                    // vectory of survey selectivty parameters on arithmetic scale

    // Fishing mortality
    init_number          log_avg_F(ph_avg_F);                    // Log average fishing mortality
    init_bounded_vector  log_F_devs(styr,endyr,-15.,15.,ph_Fdev); // Annual fishing mortality deviations
    vector               Fmort(styr,endyr);                      // Fishing mortality by year
    matrix               Z(styr,endyr,1,nages_M);                // Total mortality by year and age
    matrix               F(styr,endyr,1,nages_M);                // Fishing mortality by year and age
    matrix               S(styr,endyr,1,nages_M);                // Survivorship by year and age

    // Create a vector of natual mortalities for proj.dat
    vector               natmortv(1,nages_M);

    // Numbers at age
    init_bounded_vector  log_rec_dev(styr_rec+1,endyr,-10.,10.,ph_recdev); // Recruitment deviations from before the asssessment starts to present
    matrix               natage(styr,endyr,1,nages_M);           // Matrix of numbers at age from start year to end year

    // Biomass at age
    matrix               batage(styr,endyr,1,nages_M);           // Matrix of biomass at age from start year to end year

    // Catch at age
    matrix               catage(styr,endyr,1,nages_M)            // Matrix of predicted catch at age from start year to endyear
    vector               pred_catch_early(styr,yr_catchwt)       // Vector of predicted catches
    vector               pred_catch_later(yr_catchwt+1,endyr)    // Vector of predicted catches

    // Predicted values
    init_number          log_q_srv1(ph_q_srv1);                  // Estimate Log survey catchability
    init_number          logm(ph_m);                             // Estimate log natural mortality
    vector               pred_srv1(1,nyrs_srv1);                 // Predicted survey
    matrix               eac_fish(1,nyrs_fish_age,1,nages_D)     // Expected proportion at age in fish
    matrix               eac_srv1(1,nyrs_srv1_age,1,nages_D)     // Expected proportion at age in survey
    matrix               esc_fish(1,nyrs_fish_size,1,nlenbins)   // Expected proportion at size in fishery
    matrix               esc_srv1(1,nyrs_srv1_size,1,nlenbins)   // Expected proportion at size in survey
    
    // Effective N and SDNR    
    vector               effn_fish_age(1,nyrs_fish_age)          // Effective N for fishery age
    vector               sdnr_fish_age(1,nyrs_fish_age)          // SDNR for fishery age
    vector               effn_fish_size(1,nyrs_fish_size)        // Effective N for fishery size
    vector               sdnr_fish_size(1,nyrs_fish_size)        // SDNR for fishery size
    vector               effn_srv1_age(1,nyrs_srv1_age)          // Effective N for survey 1 age
    vector               sdnr_srv1_age(1,nyrs_srv1_age)          // SDNR for survey 1 age
    vector               effn_srv1_size(1,nyrs_srv1_size)        // Effective N for survey 1 size
    vector               sdnr_srv1_size(1,nyrs_srv1_size)        // SDNR for survey 1 size  
    
    // Standard deviation estimates for some estimated parameters
    sdreport_vector      tot_biom(styr,endyr);                   // Standard deviation report vector of total biomass
    sdreport_number      q_srv1;                                 // " " for Survey1 catchability
    sdreport_vector      pred_rec(styr,endyr);                   // " " for predicted recruitments
    vector               expl_rate(styr,endyr);                  // " " for exploitation rate 
    sdreport_number      avg_rec;                                // " " for Average recruitment 
    sdreport_number      spbiom_trend;                           // " " of Trend in biomass over last 2 years (B(t)/B(t-1); t=endyr)
    number               Depletion;                              // " " for Depletion
    sdreport_vector      spawn_biom(styr,endyr);                 // " " for spawning biomass vector
    number               natmort;                                // " " for natural mortality
    sdreport_number      LMR;
    sdreport_number      cigar;
    sdreport_number      nattymort;

    // Parameters for computing SPR rates 
    init_bounded_number  mF50(0.01,1.,ph_F50)                    // Estimated F50
    init_bounded_number  mF40(0.01,1.,ph_F50)                    // Estimated F40
    init_bounded_number  mF35(0.01,1.,ph_F50)                    // Estimated F35
    sdreport_number      F50;                                    // Standard deviation report for F50
    sdreport_number      F40;                                    // " " " F40
    sdreport_number      F35;                                    // " " " F35
    number               SB0                                     // Spawning biomass at no fishing
    number               SBF50                                   // " " at F50
    number               SBF40                                   // " " at F40
    number               SBF35                                   // " " at F35
    number               sprpen                                  // Likelihood penalty to make ADMB estimate spr rates
    matrix               Nspr(1,4,1,nages_M)                     // Matrix of number of spawners at age at each fishing mortality level

    // Likelihoods and penalty functions
    vector               surv_like(1,3);                         // Likelihood values for survey biomasses, allowance for up to 3 surveys
    vector               age_like(1,6);                          // Likelihood values for age and size compositions allowance for up 6 comps
    vector               offset(1,6);                            // Multinomial "offset"
    number               rec_like;                               // Likelihood value for recruitments
    number               ssqcatch;                               // Likelihood value for catch estimation
    number               F_mort_regularity;                      // Penalty value for fishing mortality regularity
    number               avg_sel_penalty;                        // Penalty value for selectivity regularity penalty

    // Priors
    vector               priors(1,5);                            // Prior penalty values for sigr,q,natural mortality

    // Define an objective function
    number               Like;                                   // Likelihood for data fits
    objective_function_value obj_fun;                            // Total likelihood for objective function value
    vector               pred_catch(styr,endyr);
    vector               obs_catch(styr,endyr);

    // Population projection
    matrix               N_proj(endyr+1,endyr+15,1,nages_M);
    number               FABC_proj;
    vector               FABC_tot_proj(1,nages_M);
    number               FOFL_proj;
    vector               FOFL_tot_proj(1,nages_M);
    sdreport_number      ABC;                                    // Estimate of next year's ABC
    sdreport_number      B40;
    number               OFL;
    vector               Z_proj(1,nages_M);
    vector               ZOFL_proj(1,nages_M);
    vector               S_proj(1,nages_M);
    matrix               catage_proj(endyr+1,endyr+15,1,nages_M);
    matrix               catage_proj_OFL(endyr+1,endyr+15,1,nages_M);
    vector               pred_catch_proj(endyr+1,endyr+15);
    vector               pred_catch_proj_OFL(endyr+1,endyr+15);
    sdreport_vector      spawn_biom_proj(endyr+1,endyr+15);
    sdreport_vector      tot_biom_proj(endyr+1,endyr+15);
    number               stdev_rec;
    number               FOFL;
    number               FABC;
    number               FOFL2;
    number               FABC2;
    
PROCEDURE_SECTION  

  l=l+1;                                                       // Initiate counter for random seeds in projection
  Get_Selectivity();                                           // Call function to get selectivities
  Get_Mortality_Rates();                                       // Call function to get fishing and natural mortality
  Get_Numbers_At_Age();                                        // Call function to get numbers at age per year
  Get_Catch_at_Age();                                          // Call function to get catch at age per year
  Get_Predicted_Values();                                      // Get predicted values for catch, survbio, age and size comps
  if (last_phase()){
    Get_Dependent_Vars();                                      // Solve for dependent variables like total bio, recruitment etc.
    Compute_SPR_Rates();                                       // Compute f40 etc.
    Get_Population_Projection();                               // Get 15 year population projection
  }
  Evaluate_Objective_Function();                               // Minimize objective function value
 if (mceval_phase()) {											// For outputting MCMC simulations in text format
     evalout << sigr << " " << q_srv1 << " " << F40 << " " << natmort << " " << 
     spawn_biom_proj(endyr+1) << " " <<ABC<< " " << obj_fun<< " " << tot_biom<< " " << 
     log_rec_dev<< " " << spawn_biom<< " " << log_mean_rec<< " " << spawn_biom_proj << 
     " " << pred_catch_proj<< " " << N_proj(endyr+1,1) << " " << N_proj(endyr+2,1) << " " 
     << N_proj(endyr+3,1) << " " << N_proj(endyr+4,1) << " " << N_proj(endyr+5,1) << " " 
     << N_proj(endyr+6,1) << " " << N_proj(endyr+7,1) << " " << N_proj(endyr+8,1) << " " 
     << N_proj(endyr+9,1) << " " << N_proj(endyr+10,1) << " " << tot_biom_proj << "\n";
  }

// flag - has hard coded recruitment age (4) starting at j=1 and using double(j+3)
FUNCTION Get_Selectivity
    for (j=1;j<=nages_M;j++){
      fish_sel(j) = 1./(1. + mfexp(-2.944438979 * (double(j+3) - a50) / delta));
      srv1_sel(j) = 1./(1. + mfexp(-2.944438979 * (double(j+3) - a50_srv1) / delta_srv1));
    }

FUNCTION Get_Mortality_Rates
    natmort = mfexp(logm);			// setting natural mortality to arithmetic scale
    if(ph_m>0) nattymort = natmort;
       else nattymort = log_mean_rec;
    Fmort = mfexp(log_avg_F + log_F_devs);	//setting fishing mortaltiy to arithmetic scale
    for (iyr=styr; iyr<=endyr; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel;		// Getting fully selected fishing mortality
    Z = F + natmort;				// Fully selected total mortality
    S = mfexp(-1.0 * Z);

FUNCTION Get_Numbers_At_Age  
    // Start year abundance
    int itmp;
    for (j=2;j<nages_M;j++) {
        itmp = styr+1-j;
        natage(styr,j) = mfexp(log_mean_rec - natmort * double(j-1)+ log_rec_dev(itmp)); 
    }
    natage(styr,nages_M) = mfexp(log_mean_rec - natmort * (nages_M-1)) / (1. - exp(-natmort));

    // Remaining years
    for (i=styr;i<endyr;i++) {
        natage(i,1) = mfexp(log_rec_dev(i) + log_mean_rec);
        natage(i+1)(2,nages_M) = ++elem_prod(natage(i)(1,nages_M-1),S(i)(1,nages_M-1));       // Following year
        natage(i+1,nages_M) += natage(i,nages_M) * S(i,nages_M);
        Sp_Biom(i) = natage(i) * wt_mature;                        // Old way, correct way to get SSB at start of spawning season would be: Sp_Biom(i) = elem_prod(natage(i),pow(S(i),spawn_fract)) * wt_mature;
    }

    // End year abundance
    natage(endyr,1) = mfexp(log_rec_dev(endyr) + log_mean_rec); 
    // BETTER WAY Sp_Biom(endyr) = elem_prod(natage(endyr),pow(S(endyr),spawn_fract)) * wt_mature;  //Right way, old way was: Sp_Biom(endyr) = natage(endyr)* wt_mature;
    Sp_Biom(endyr) = natage(endyr)* wt_mature; //Old way

FUNCTION Get_Catch_at_Age
    pred_catch_early.initialize();
    pred_catch_later.initialize();
    for (iyr=styr; iyr<=yr_catchwt; iyr++){
        catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
        pred_catch_early(iyr) = catage(iyr) * wt;
    }
    for (iyr=yr_catchwt+1; iyr<=endyr; iyr++) {
        catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
        pred_catch_later(iyr) = catage(iyr) * wt;
    }

FUNCTION Get_Predicted_Values    
    offset.initialize();

    // Calculate predicted data values
    q_srv1 = exp(log_q_srv1); // Survey catchability at arithmetic scale
    for (i=1;i<=nyrs_srv1;i++){
        pred_srv1(i) = q_srv1 * (natage(yrs_srv1(i))*elem_prod(srv1_sel,wt));   // Predicted Survey biomass
        }

    // Predicted Fishery age comps, N, effn, sdnr, and offset 
    for (i=1;i<=nyrs_fish_age;i++) {
    eac_fish(i)  = catage(yrs_fish_age(i))/sum(catage(yrs_fish_age(i))) * ageage(age_age_ind_fsh(i));
    //nmulti_fish_age(i) = sqrt(nsamples_fish_age(i));
    effn_fish_age(i) = (1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i));
    sdnr_fish_age(i) = sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)));
    }
    for (i=1; i<=nyrs_fish_age; i++) {
    oac_fish(i)/=sum(oac_fish(i));
    offset(1) -= nmulti_fish_age(i) *((oac_fish(i) + 0.00001)*log(oac_fish(i) + 0.00001)); 
    }

    // Predicted Survey1 age comps, N, effn, sdnr, and offset       
    for (i=1;i<=nyrs_srv1_age;i++) {
    eac_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_age(i)))/(natage(yrs_srv1_age(i)) * srv1_sel)* ageage(age_age_ind_srv(i));
    //nmulti_srv1_age(i) = sqrt(nsamples_srv1_age(i));
    effn_srv1_age(i) = (1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i));
    sdnr_srv1_age(i) = sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)));
    }
    for (i=1; i<=nyrs_srv1_age; i++) {
    oac_srv1(i)/=sum(oac_srv1(i));
    offset(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001)*log(oac_srv1(i) + 0.00001));
    }

    // Predicted Fishery size comps, N, effn, sdnr, and offset
    for (i=1;i<=nyrs_fish_size;i++) {
    esc_fish(i)  = catage(yrs_fish_size(i))/sum(catage(yrs_fish_size(i))) * sizeage(siz_age_ind_fsh(i));
    //nmulti_fish_size(i) = nhauls_fish_size(i);
    effn_fish_size(i) = (1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i));
    sdnr_fish_size(i) = sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)));
    }
    //nmulti_fish_size = nmulti_fish_size/max(nmulti_fish_size)*100;
    for (i=1; i<=nyrs_fish_size; i++) {
    osc_fish(i)/=sum(osc_fish(i));
    offset(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001)*log(osc_fish(i) + 0.00001));
    }

    // Predicted Survey1 size comps, N, effn, sdnr, and offset
    for (i=1;i<=nyrs_srv1_size;i++) {
    esc_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_size(i))) /(natage(yrs_srv1_size(i)) * srv1_sel)* sizeage(siz_age_ind_srv1(i));
    effn_srv1_size(i) = (1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i));
    //sdnr_srv1_size(i) = sdnr(esc_srv1(i),osc_srv1(i),double(effn_srv1_size(i)));
    }
    for (i=1; i<=nyrs_srv1_size; i++) {
    osc_srv1(i)/=sum(osc_srv1(i));
    offset(4) -= effn_srv1_size(i)*((osc_srv1(i) + 0.00001)*log(osc_srv1(i) + 0.00001));
    }

    // Predicted catch
    pred_catch(styr,yr_catchwt) = pred_catch_early;
    pred_catch(yr_catchwt+1,endyr) = pred_catch_later;
    obs_catch(styr,yr_catchwt) = obs_catch_early;
    obs_catch(yr_catchwt+1,endyr) = obs_catch_later;
    
    // set up some sdreport numbers
    cigar = sigr;
    LMR = log_mean_rec;

FUNCTION Get_Dependent_Vars

  for (i=styr;i<=endyr;i++) {
    pred_rec(i) = natage(i,1);                                 // Setting up results based on estimated paramters
    tot_biom(i) = wt * natage(i);                              // Total biomass results
    expl_rate(i) = pred_catch(i)/tot_biom(i);                  // Setting up results based on estimated paramters
    spawn_biom(i) = Sp_Biom(i) ;                               // Spawning biomass result
  }
  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);         // 1-Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);

FUNCTION Compute_SPR_Rates
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0 = 0.;
  SBF50 = 0.;
  SBF40 = 0.;
  SBF35 = 0.;
  // Scale F-spr rates to be on full-selected values
  F50 = mF50 * max(fish_sel);
  F40 = mF40 * max(fish_sel);
  F35 = mF35 * max(fish_sel);
  for (i=1;i<=4;i++){
    Nspr(i,1)=1.;
  }
  for (j=2;j<nages_M;j++) {
    Nspr(1,j) = Nspr(1,j-1) * mfexp(-1. * natmort);
    Nspr(2,j) = Nspr(2,j-1) * mfexp(-1. * (natmort + mF50 * fish_sel(j-1)));
    Nspr(3,j) = Nspr(3,j-1) * mfexp(-1. * (natmort + mF40 * fish_sel(j-1)));
    Nspr(4,j) = Nspr(4,j-1) * mfexp(-1. * (natmort + mF35 * fish_sel(j-1)));
  }
  Nspr(1,nages_M) = Nspr(1,nages_M-1) * mfexp(-1. * natmort) / (1.-mfexp(-1.*natmort));
  Nspr(2,nages_M) = Nspr(2,nages_M-1) * mfexp(-1. * (natmort + mF50 * fish_sel(nages_M-1))) /
                        (1.-mfexp(-1. * (natmort + mF50 * fish_sel(nages_M))));
  Nspr(3,nages_M) = Nspr(3,nages_M-1) * mfexp(-1. * (natmort + mF40 * fish_sel(nages_M-1))) /
                        (1.-mfexp(-1. * (natmort + mF40 * fish_sel(nages_M))));
  Nspr(4,nages_M) = Nspr(4,nages_M-1) * mfexp(-1. * (natmort + mF35 * fish_sel(nages_M-1))) /
                        (1.-mfexp(-1. * (natmort + mF35 * fish_sel(nages_M))));
  
  for (j=1; j<=nages_M; j++) {
   // Kill them off till (spawn_fract)
    SB0    += Nspr(1,j) * wt_mature(j) * mfexp(-spawn_fract * natmort);
    SBF50  += Nspr(2,j) * wt_mature(j) * mfexp(-spawn_fract * (natmort + mF50 * fish_sel(j)));
    SBF40  += Nspr(3,j) * wt_mature(j) * mfexp(-spawn_fract * (natmort + mF40 * fish_sel(j)));
    SBF35  += Nspr(4,j) * wt_mature(j) * mfexp(-spawn_fract * (natmort + mF35 * fish_sel(j)));
   }
  sprpen    = 100. * square(SBF50 / SB0-0.5);
  sprpen   += 100. * square(SBF40 / SB0-0.4);
  sprpen   += 100. * square(SBF35 / SB0-0.35);
  B40 = SBF40 * mean(pred_rec(1979,endyr-recage));

FUNCTION Get_Population_Projection
  //  Abundance at start of first projection year
  int k;
// Flag ----------------------------------
  // Recruitment in endyr+1
//  if(mceval_phase()) {
//    stdev_rec = sqrt(norm2(value(log_rec_dev(1979 + recage,endyr - recage))-mean(value(log_rec_dev(1979 + recage,endyr - recage))))/(size_count(value(log_rec_dev(1979 + recage,endyr - recage))) - 1));
//    k=round(value(stdev_rec) * 10000);
//    N_proj(endyr+1,1) = mfexp(value(log(mean(value(pred_rec(1979+recage,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));
//    cout << stdev_rec<< " " <<k<< " " << l<< " " << "\n";
//  }  else {
//    N_proj(endyr+1,1)= value(mean(pred_rec(1979+recage,endyr-recage))); 
//  }

  if(mceval_phase()) {
    stdev_rec = sqrt(norm2(value(log_rec_dev(1979, endyr - recage))-mean(value(log_rec_dev(1979, endyr - recage))))/(size_count(value(log_rec_dev(1979, endyr - recage))) - 1));
    k=round(value(stdev_rec) * 10000);
    N_proj(endyr+1,1) = mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));
    cout << stdev_rec<< " " <<k<< " " << l<< " " << "\n";
  } else {
    N_proj(endyr+1,1)= value(mean(pred_rec(1979, endyr-recage))); 
    }

  // Abundance for remaining age classes in endyr+1
  for (j=1;j<nages_M-1;j++)
    N_proj(endyr+1,j+1) = natage(endyr,j) * S(endyr,j);
  N_proj(endyr+1,nages_M) = natage(endyr,nages_M-1) * S(endyr,nages_M-1) + natage(endyr,nages_M) * S(endyr,nages_M);
  tot_biom_proj(endyr+1) = N_proj(endyr+1) * wt;
  spawn_biom_proj(endyr+1) = elem_prod(N_proj(endyr+1),pow(mfexp(-yieldratio * FABC_tot_proj-natmort),spawn_fract)) * wt_mature;

  // loop through to endyr+15
  for (i=endyr+1;i<=endyr+15;i++) {

   // F ABC 
    if (spawn_biom_proj(i)/B40 > 1.) {
      FABC_proj = F40;
      FOFL_proj = F35; 
      } else {
      FABC_proj = F40 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); 
      FOFL_proj = F35 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); 
      }
    for (j=1;j<=nages_M;j++) {  
      FOFL_tot_proj(j) = fish_sel(j) * FOFL_proj;
      FABC_tot_proj(j) = fish_sel(j) * FABC_proj;
      Z_proj(j) = FABC_tot_proj(j) + natmort;
      ZOFL_proj(j) = FOFL_tot_proj(j) + natmort;
      S_proj(j) = mfexp(-1.0 * Z_proj(j));
    }

   // Catch 
    for (j=1;j<=nages_M;j++) { 
      catage_proj(i,j) = yieldratio * N_proj(i,j) * FABC_tot_proj(j) / Z_proj(j) * (1.-S_proj(j));
      catage_proj_OFL(i,j) = yieldratio * N_proj(i,j) * FOFL_tot_proj(j) / ZOFL_proj(j) * (1. - mfexp(-ZOFL_proj(j)));
    }
    pred_catch_proj(i) = catage_proj(i) * wt / yieldratio;
    pred_catch_proj_OFL(i) = catage_proj_OFL(i) * wt / yieldratio;
   
   // Next year's abundance
    if (i < endyr + 15) {
     if (mceval_phase()) {
       stdev_rec = sqrt(norm2(value(log_rec_dev(1979+recage,endyr-recage))-mean(value(log_rec_dev(1979+recage,endyr-recage))))/(size_count(value(log_rec_dev(1979+recage,endyr-recage)))-1));
       k=round(value(spawn_biom(endyr)*10000))+i;
       k=k+i;
       N_proj(i+1,1)= mfexp((log(mean(value(pred_rec(1979+recage,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l))); }
     else {
       N_proj(i+1,1)= value(mean(pred_rec(1979,endyr-recage))); }
     for (j=1; j<nages_M-1;j++) {
       N_proj(i+1,j+1) = N_proj(i,j) * mfexp(-yieldratio*FABC_tot_proj(j)-natmort); 
       }
     N_proj(i+1,nages_M) = N_proj(i,nages_M-1) * mfexp(-yieldratio*FABC_tot_proj(nages_M-1)-natmort)+ N_proj(i,nages_M) * mfexp(-yieldratio*FABC_tot_proj(nages_M)-natmort);
     spawn_biom_proj(i+1) = elem_prod(N_proj(i+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;
     tot_biom_proj(i+1) = N_proj(i+1)*wt;
     }
   }

  // Set up control rules for harvest rates
  if (spawn_biom_proj(endyr+1)/B40 > 1.) {
      FABC = F40;
      FOFL = F35; 
      FABC2 = F40;
      FOFL2 = F35;
  } else {
      FABC = F40 * (spawn_biom_proj(endyr+1) / B40 - 0.05) / (1 - 0.05); 
      FOFL = F35 * (spawn_biom_proj(endyr+1) / B40 - 0.05) / (1 - 0.05);  
      FABC2 = F40 * (spawn_biom_proj(endyr+2) / B40 - 0.05) / (1 - 0.05); 
      FOFL2 = F35 * (spawn_biom_proj(endyr+2) / B40 - 0.05) / (1 - 0.05);
  }

  // Projected ABC and OFL
  OFL = pred_catch_proj_OFL(endyr+1);
  ABC = pred_catch_proj(endyr+1);

FUNCTION Evaluate_Objective_Function   
  Like.initialize();
  // Call functions to compute data likelihood
  Catch_Like();                                                // Catch biomass likelihood (lognormal)
  Surv_Like();                                                 // Trawl survey biomass likelihood (lognormal)
  Size_Age_Like();                                             // Age/Size composition likelihood (multinomial)
  
  // Call functions to compute prior penalties
  Calc_Priors();                                               // Prior penalties for estimated parameters
  Rec_Like();                                                  // Penalty function for recruitment
  F_Like();                                                    // Penalty function for fishing mortality deviations

  // Sum objective function
  Like              += ssqcatch ;
  Like              += sum(surv_like);
  Like              += sum(age_like);
  obj_fun           += Like;
  obj_fun           += sum(priors);
  obj_fun           += wt_rec_var * rec_like;
  if(active(log_F_devs))
    obj_fun         += F_mort_regularity;
  if (current_phase()<3)
      obj_fun       += norm2(F);   
  if (active(mF50)&&last_phase())
    obj_fun         += sprpen;                                 // To solve for the F40 etc.     

FUNCTION Catch_Like
  ssqcatch.initialize();
  ssqcatch  +=  wt_ssqcatch *norm2(log(obs_catch_early+.00001)-log(pred_catch_early+.00001));
  ssqcatch  +=  wt_ssqcatch2 *norm2(log(obs_catch_later+.00001)-log(pred_catch_later+.00001));

FUNCTION Surv_Like
  surv_like.initialize();
  // Calculate likelihood for survey biomass
  for (i=1; i<=nyrs_srv1; i++) {
	  	  surv_like(1) += square(obs_srv1_biom(i)-pred_srv1(i) )/ (2.*square(obs_srv1_se(i))); }  // likelihood for survey biomass

  surv_like(1) *= wt_srv1 ;  

FUNCTION Size_Age_Like
  age_like.initialize();
  // Calculate multinomial likelihoods for survey age, fishery size, and survey size and subtract "offset"
  for (i=1; i <= nyrs_fish_age; i++)
    age_like(1) -= nmulti_fish_age(i) * ((oac_fish(i) + 0.00001) * log(eac_fish(i) + 0.00001));
  age_like(1)   -= offset(1);

  for (i=1; i <= nyrs_srv1_age; i++)
    age_like(2) -= nmulti_srv1_age(i) * ((oac_srv1(i) + 0.00001) * log(eac_srv1(i) + 0.00001));
  age_like(2)   -= offset(2);

  for (i=1; i <= nyrs_fish_size; i++)
    age_like(3) -= nmulti_fish_size(i) * ((osc_fish(i) + 0.00001) * log(esc_fish(i) + 0.00001));
  age_like(3)   -= offset(3);

  for (i=1; i <= nyrs_srv1_size; i++)
    age_like(4) -= effn_srv1_size(i) * ((osc_srv1(i) + 0.00001) * log(esc_srv1(i) + 0.00001));
  age_like(4)   -= offset(4);

 
  // Multiple each likelihood by their weights from .ctl file
  age_like(1) *= wt_fish_age;
  age_like(2) *= wt_srv1_age;
  age_like(3) *= wt_fish_size;
  age_like(4) *= wt_srv1_size;

FUNCTION Calc_Priors
  // Calculate prior penalties
    priors.initialize();
    if (active(sigr))
      priors(1)    = square(log(sigr/sigrprior)) / (2.*square(cvsigrprior));
    if (active(log_q_srv1))
      priors(2)    = square(log_q_srv1-log_q_srv1prior) / (2.*square(cvq_srv1prior));
    if (active(logm))
      priors(3)    = square(logm-log(mprior)) / (2.*square(cvmprior));
FUNCTION Rec_Like
  rec_like.initialize();

      //if (active(sigr)){
        //rec_like = norm2(log_rec_dev+sigr*sigr/2.) / (2.*square(sigr)) + size_count(log_rec_dev) * log(sigr);
      //} else { 
        rec_like = norm2(log_rec_dev)/(2*square(sigr)) + (size_count(log_rec_dev)*log(sigr));
        //rec_like = norm2(log_rec_dev+sigr*sigr/2.) / (2.*square(sigr)) ;
      //}

FUNCTION F_Like
   
  F_mort_regularity.initialize();

  if(active(log_F_devs))
    F_mort_regularity  = wt_fmort_reg * norm2(log_F_devs);

FUNCTION double round(double r) 
    return double((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5)); 
REPORT_SECTION
  cout << "-------------Finished: " << current_phase() << " " << Like << " " << age_like << "\n";

  // Beginning of all outputting
  // Goes to routine that automatically creates input file for projection model
  if (last_phase())
    write_proj();

  // Output file (.rep) which is loaded into R for data outputs
  report << "~~~~~~~ Executive Summary Material ~~~~~~~" << "\n";report << "\n";
  report << "     Model name" << "\n";
  report << model_name << "\n";
  report << "     .dat file" << "\n";
  report << data_file << "\n";
  report << "     Number parameters estimated" << "\n";
  report << initial_params::nvarcalc() << "\n";
  report << "     TotalBiomass for " <<endyr+1<< "\n";
  report << N_proj(endyr+1)*wt << "\n";
  report << "     TotalBiomass for " <<endyr+2<< "\n";
  report << N_proj(endyr+2)*wt << "\n";
  report << "     Female_Spawning Biomass for " <<endyr+1<< "\n";
  report << spawn_biom_proj(endyr+1) << "\n";
  report << "     Female_Spawning_Biomass for " <<endyr+2<< "\n";
  report << spawn_biom_proj(endyr+2) << "\n";
  report << "     B_zero"     << "\n";
  report << SB0 * mean(pred_rec(1977+recage,endyr-recage)) << "\n";
  report << "     B_40" << "\n";
  report << B40<< "\n";
  report << "     B_35" << "\n";
  report << SBF35 * mean(pred_rec(1977+recage,endyr-recage)) << "\n";
  report << "     F_40" << "\n";
  report << F40<< "\n";
  report << "     F_35" << "\n";
  report << F35<< "\n";
  report << "     F_ABC for " <<endyr+1 << "\n";
  report << FABC<< "\n";
  report << "     F_ABC for " <<endyr+2 << "\n";
  report << FABC2<< "\n";
  report << "     ABC for " <<endyr+1 << "\n";
  report << pred_catch_proj(endyr+1) << "\n";
  report << "     ABC for " <<endyr+2 << "\n";
  report << pred_catch_proj(endyr+2) << "\n";
  report << "     F_OFL for " <<endyr+1 << "\n";
  report << FOFL << "\n";
  report << "     F_OFL for " <<endyr+2 << "\n";
  report << FOFL2 << "\n";
  report << "     OFL for " <<endyr+1 << "\n";
  report << OFL << "\n"; 
  report << "     OFL for " <<endyr+2 << "\n";
  report << pred_catch_proj_OFL(endyr+2) << "\n"; 
  report << "     Total likelihood" << "\n";
  report << obj_fun << "\n";
  report << "     Data likelihood" << "\n";
  report << Like << "\n" << "\n";
  
  report << "~~~~~~~ Some more key parameter estimates ~~~~~~~" << "\n" << "\n";
  report << "   q_trawl" << "\n";
  report << q_srv1<< "\n";
  report << "   nat_mort" << "\n";
  report << natmort << "\n";
  report << "  sigr" << "\n";  
  report << sigr<< "\n";  
  report << "   log_mean_rec" << "\n";
  report << log_mean_rec<< "\n";

  report << "~~~~~~~ Rest of model/data output ~~~~~~~" << "\n" << "\n";
  report << "Year " << yy << "\n";
  report << "Pred_Catch " << pred_catch_early << pred_catch_later << "\n";
  report << "Obs_Catch " << obs_catch_early<< obs_catch_later << "\n" << "\n";
  report << "Catch_at_age " << aa << "\n";
  for (i=styr;i<=endyr;i++) report << i << " " << catage(i) << "\n";
  report << "\n";
  report << "Numbers " << aa << "\n";
  for (i=styr;i<=endyr;i++) report << i << " " << natage(i) << "\n"; 
  report << "\n";
  
  report << "Obs_P_fish_age" << aa_D << "\n";
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i) << " " << oac_fish(i) << "\n";
  report << "\n";
  report << "Pred_P_fish_age" << aa_D << "\n";
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i) << " " << eac_fish(i) << "\n";
  report << "\n";
  report << "#-- yrs_fish_age nmulti_fish_age effn_fish_age sdnr_fish_age" << "\n";
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i) << " " << nmulti_fish_age(i) << " " << effn_fish_age(i) << " " << sdnr_fish_age(i)  << " #!!" << "\n";
  report << "\n";
  report << "Obs_P_fish_size" << len_bin_labels << "\n";
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i) << " " << osc_fish(i) << "\n";
   report << "\n";
  report << "Pred_P_fish_size" << len_bin_labels << "\n";
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i) << " " <<esc_fish(i) << "\n";
  report << "\n";
  report << "#-- yrs_fish_size nmulti_fish_size effn_fish_size sdnr_fish_size" << "\n";
  for (i=1;i<=nyrs_fish_size;i++) report <<yrs_fish_size(i) << " " <<nmulti_fish_size(i) << " " <<effn_fish_size(i) << " " << sdnr_fish_size(i) << " #!!" << "\n";
  
  report << "\n";
   report << "Obs_P_srv1_age" << aa_D << "\n";
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i) << " " << oac_srv1(i) << "\n";
  report << "\n";  
  report << "Pred_P_srv1_age" << aa_D << "\n";
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i) << " " <<eac_srv1(i) << "\n";
  report << "\n";  
  report << "#-- yrs_srv1_age nmulti_srv1_age effn_srv1_age sdnr_srv1_age" << "\n";
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i) << " " << nmulti_srv1_age(i) << " " <<effn_srv1_age(i) << " " << sdnr_srv1_age(i) << " #!!" << "\n";
  report << "\n";  
  report << "Obs_P_srv1_size" << len_bin_labels << "\n";
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i) << " " << osc_srv1(i) << "\n";
  report << "\n";  
  report << "Pred_P_srv1_size" << len_bin_labels << "\n";
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i) << " " <<esc_srv1(i) << "\n";
  report << "\n";  
  report << "Bottom Trawl Survey Biomass " << "\n";
  report << "Year:     " << yrs_srv1 << "\n";
  report << "Predicted:   " << pred_srv1 << "\n";
  report << "Observed:   " << obs_srv1_biom << "\n";
  report << "Observed_SE:   " << obs_srv1_se << "\n" << endl;
  report << "Year " << yy << endl;
  report << "SpBiom " << spawn_biom << "\n";
  report << "Tot_biom " << tot_biom   << "\n";
  report << "Recruitment " << pred_rec   << "\n";
  report << "Fully_selected_F " << Fmort * max(fish_sel) << "\n" << "\n";;
  report << "Age  " << aa<< endl;
  report << "Weight " << wt << endl;
  report << "Maturity " << p_mature << "\n";
  report << "Fishery_Selectivity_1967-1976" << fish_sel  << "\n";
  report << "Bottom_Trawl_Survey_Selectivity " << srv1_sel / max(srv1_sel) << "\n";
    report << "F35 F40 F50 " << "\n";
  report <<  F35 << " " << F40 << " " <<  F50 << "\n";report << "\n";
  report << "~~~~~~~ Wts and Likelihoods ~~~~~~~" << endl;report << "\n";
  report << wt_ssqcatch << " " << ssqcatch  << " SSQ Catch Likelihood" << endl;
  report << wt_srv1 << " " << surv_like(1) << " Bottom Trawl Survey Likelihood" << endl;
  report << wt_fish_age << " " << age_like(1)  << " Fishery Age Composition Likelihood"  << endl;
  report << wt_srv1_age << " " << age_like(2)  << " Bottom Trawl Survey Age Composition Likelihood" << endl;
  report << wt_fish_size << " " << age_like(3)  << " Fishery Size Composition Likelihood" << endl;
  report << wt_srv1_size << " " << age_like(4)  << " Bottom Trawl Survey Size Composition Likelihood" << endl;
  report << wt_rec_var << " " << rec_like     << " Recruitment Deviations Likelihood" << endl;
  report << wt_fmort_reg << " " << F_mort_regularity<< " Fishing Mortality Deviations Penalty" << endl;
  report << 1 << " " << priors(1)  << " Priors SigmaR" << "\n";
  report << 1 << " " << priors(2)  << " Priors q Bottom Trawl Survey" << "\n";
  report << 1 << " " << priors(3)  << " Priors M" << "\n";
  report << 1 << " " << sprpen  << "  Spawner-recruit penalty" << "\n";
  report << " " << obj_fun    << "  Objective Function" << "\n";
  report << " " << Like       << "  Data Likelihood" << "\n" << "\n";
  report << "~~~~~~~ Projection outputs ~~~~~~~" << "\n" << "\n";
  report << "N_at_age projected " << "\n" 
          << N_proj << "\n" 
          << " spawn_bio projected" << "\n" 
          << spawn_biom_proj << "\n";

FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;

FUNCTION write_proj
 ofstream newproj("proj.dat");
// Function to write out data file for new Ianelli 2005 projection model....
 newproj << "#Species name here:" << "\n";
 newproj << model_name+"_"+data_file << "\n";
 newproj << "#SSL Species?" << "\n";
 newproj << "0" << "\n";
 newproj << "#Constant buffer of Dorn?" << "\n";
 newproj << "0" << "\n";
 newproj << "#Number of fisheries?" << "\n";
 newproj << "1" << "\n";
 newproj << "#Number of sexes?" << "\n";
 newproj << "1" << "\n";
 newproj << "#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)" << "\n";
 newproj << mean(Fmort(endyr-4,endyr)) << "\n";
 newproj << "#_Author_F_as_fraction_F_40%" << "\n";
 newproj << "1" << "\n";
 newproj << "#ABC SPR" << "\n";
 newproj << "0.4" << "\n";
 newproj << "#MSY SPR" << "\n";
 newproj << "0.35" << "\n";
 newproj << "#_Spawn_month" << "\n";
 newproj << spawn_fract*12+1<< "\n";
 newproj << "#_Number_of_ages" << "\n";
 newproj <<nages_M<< "\n";
 newproj << "#_F_ratio(must_sum_to_one_only_one_fishery)" << "\n";
 newproj << "1" << "\n";
 for (j=1;j<=nages_M;j++) natmortv = natmort; 
 newproj << "#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<< "\n";
 newproj << "#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2" << aa<< "\n" << p_mature << endl;
 newproj << "#_Wt_at_age_spawners" << aa<< "\n" <<wt << endl;
 newproj << "#_Wt_at_age_fishery" << aa<< "\n" <<wt << endl;
 newproj << "#_Selectivity_fishery_scaled_to_max_at_one" << aa<< "\n" <<fish_sel/max(fish_sel) << endl;
 newproj << "#_Numbers_at_age_end_year" << aa<< "\n" <<natage(endyr) << endl;
 newproj << "#_N_recruitment_years" << "\n" <<endyr-recage-1979+1<< endl;
 newproj << "#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits" <<yy(1979,endyr-recage) << "\n" << pred_rec(1979,endyr-recage) << endl;
 newproj << "#_Spawners per recruitment (starting at 1977)" << "\n" << spawn_biom(1977,endyr-recage)/1000<< endl;
 newproj.close();

GLOBALS_SECTION

 # include "admodel.h"
  adstring model_name;
  adstring data_file;
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
