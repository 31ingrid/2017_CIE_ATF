//                      Kam.tpl
//           Uses the shelf survey, slope survey and the Aleutian Islands (fits biomass and length and age comps)
//           domed shaped selectivity for the shelf survey males,
//           can estimate the temperature effect on shelf survey catchability
DATA_SECTION
!!CLASS ofstream evalout("kam.mcmc.out");
  init_int styr         //start year of model
  init_int endyr        //end year
  init_int styr_fut     //start year of projections (endyr+1) 
  init_int endyr_fut    //end year of projections
  init_number srv2_qinit //initial value for slope survey q
  init_number M_init    //initial value for M
  init_number phase_M  
//  init_number mfavg     //female natural mortality rate
//  init_number mfcv      //cv of female natural mortality rate
//  init_number mmavg     //male natural mortality rate
//  init_number mmcv      //cv of male natural mortality rate
//  init_number q1avg     //avg catchability
//  init_number q1cv      //cv or avg catchability
  init_int phase_F40      //phase F40 is estimated
  init_int phase_Temp      //phase Temp is estimated
  init_number median_rec  //median recruit value to use for the last 3 years
//  init_int nages_read     //# of ages read
  init_int nages          //# of ages in the model
 //selectivity is set to the selectivity at nselages-1 after age nselages 
  init_int nselages       //fishery (for asymptotic selectivity)
  init_int nselages_srv1  //slope survey (for asymptotic selectivity) 
  init_int nselages_srv2  //slope survey (for asymptotic selectivity)
  init_int nselages_srv3  //Aleutian Islands survey (for asymptotic selectivity)
//  init_int monot_sel    //selectivity smoothing function for fishery
//  init_int monot_sel_srv1  //selectivity smoothing function for survey
  init_int phase_logistic_sel
//  init_int phase_selcoffs
//  init_vector wt_like(1,8)  //used for smoothing selectivity function
//  !! phase_logistic_sel = -2;
//  !! if (phase_logistic_sel<0) phase_selcoffs=2; else phase_selcoffs=-2;
 //sample size for length comps for weighting likelihoods  
  init_int nlen             //# of length bins
  !! cout << nlen   <<endl;
  init_int nobs_fish_length          //# of years of fishery data
  init_ivector yrs_fish_length(1,nobs_fish_length)   //years with fishery data
  init_matrix nsamples_fish(1,2,1,nobs_fish_length)  //sample size (weights) for each sex and yr of fishery data
  init_int nobs_srv1          //# of years of shelf survey data
  !! cout << nobs_srv1  << endl;
  init_int nobs_srv2          //# of years of slope survey data
  !! cout <<nobs_srv2<<endl;
  init_int nobs_srv3          //# of years of Aleutian Islands data
  !! cout <<nobs_srv3<< endl; 
  init_ivector yrs_srv1(1,nobs_srv1)   //years with shelf survey data
  !! cout <<yrs_srv1<<endl;
  init_ivector yrs_srv2(1,nobs_srv2)   //years with slope survey data
  !! cout <<yrs_srv2<<endl;  
  init_ivector yrs_srv3(1,nobs_srv3)   //years with Aleutian Islands survey data
  !! cout <<yrs_srv3<<endl;
  init_int nobs_srv1_length          //# yrs with shelf survey length data
  !! cout <<nobs_srv1_length<< endl;
  init_int nobs_srv2_length          //# yrs with slope survey length data
  !! cout <<nobs_srv2_length<<endl;
  init_int nobs_srv3_length          //# yrs with Aleutian Islands survey length data
  !! cout<<nobs_srv3_length<<endl; 
  init_ivector yrs_srv1_length(1,nobs_srv1_length)    //yrs with shelf survey length data
  init_ivector yrs_srv2_length(1,nobs_srv2_length)    //yrs with slope survey length data
  init_ivector yrs_srv3_length(1,nobs_srv3_length)    //yrs with Aleutian Islands survey length data
  init_matrix nsamples_srv1_length(1,2,1,nobs_srv1_length)  //sample size for each length comp by sex and year from shelf survey
  init_matrix nsamples_srv2_length(1,2,1,nobs_srv2_length)  //sample size for each length comp by sex and year from slope survey
  init_matrix nsamples_srv3_length(1,2,1,nobs_srv3_length)  //sample size for each length comp by sex and year from Aleutian I survey
  init_3darray obs_p_srv1_length(1,2,1,nobs_srv1_length,1,nlen) //shelf survey length comps by bin, sex and yr
  init_3darray obs_p_srv2_length(1,2,1,nobs_srv2_length,1,nlen) //slope survey length comps by bin, sex and yr
  init_3darray obs_p_srv3_length(1,2,1,nobs_srv3_length,1,nlen) //Aleutian Islands survey length comps by bin, sex and yr
  !! cout << obs_p_srv1_length(1,1,1)<<endl;
  //!! cout << obs_p_srv2_length(1,1,7)<<endl;
  !! cout << obs_p_srv3_length<<endl;
  init_3darray obs_p_fish(1,2,1,nobs_fish_length,1,nlen)  //fishery length comps
  !! cout<<obs_p_fish<<endl;
  init_vector catch_bio(styr,endyr)    //catch by year
  !! cout<<  catch_bio        <<endl;
  init_vector obs_srv1(1,nobs_srv1)    //shelf survey biomass by year
  init_vector obs_srv1_sd(1,nobs_srv1) //shelf survey SE by year
  init_vector obs_srv2(1,nobs_srv2)    //slope survey biomass by year
  init_vector obs_srv2_sd(1,nobs_srv2) //slope survey SE by year
  init_vector obs_srv3(1,nobs_srv3)    //Aleutian Islands survey biomass by year
  init_vector obs_srv3_sd(1,nobs_srv3) //Aleutian Islands survey SE by year
 //need wt vector by length for split sex?
 //init_vector wt(1,nlen)        
  init_matrix wt(1,2,1,nages)          //weight-at-age by sex
  init_vector maturity(1,nages)        //female prop. mature-at-age

  //!! cout<<  maturity(1,nages)        <<endl;
  //!! cout<<  sel_srv2<<endl;

 //this is how to change the name of the data file
 //can put at the beginning of the file to read everything or
 //part way down
 // LOCAL_CALCS
 //  ad_comm::change_datafile_name("atf.ctl")
 // END_CALCS
 //Local_calcs is indented 1 space
 // LOCAL_CALCS
 //  cout<<maturity<<endl;
 // END_CALCS
//length age transition matrix

  init_3darray lenage(1,2,1,nages,1,nlen)  //length-age transition matrix
  !! cout<<  lenage(2,nages)<<endl;
  !! cout <<nobs_srv1<<endl;
  init_vector bottom_temps(1,nobs_srv1)    //shelf survey bottom temperatures
  !! cout<<bottom_temps(1,nobs_srv1)<<endl;
 init_int nobs_srv2_age                   // # of years with slope survey ages
  !! cout <<nobs_srv2_age<<endl;
  init_ivector yrs_srv2_age(1,nobs_srv2_age)  //years of slope survey with ages
  !! cout <<yrs_srv2_age(1,nobs_srv2_age)<<endl;
  init_matrix nsamples_srv2_age(1,2,1,nobs_srv2_age)   //sample size of slope ages read in each year, by sex
  !! cout <<nsamples_srv2_age<<endl;
  init_3darray obs_p_srv2_age(1,2,1,nobs_srv2_age,1,nages)  //slope survey age comps by sex and year
  !! cout <<obs_p_srv2_age<<endl; 
 init_int nobs_srv3_age                   // # of years with Aleutian survey ages
  !! cout <<nobs_srv3_age<<endl;
  init_ivector yrs_srv3_age(1,nobs_srv3_age)  //years of Aleutian  survey with ages
  !! cout <<yrs_srv3_age(1,nobs_srv3_age)<<endl;
  init_matrix nsamples_srv3_age(1,2,1,nobs_srv3_age)   //sample size of ages read in each year, by sex
  !! cout <<nsamples_srv3_age<<endl;
  init_3darray obs_p_srv3_age(1,2,1,nobs_srv3_age,1,nages)  //Aleutian survey age comps by sex and year
  !! cout <<obs_p_srv3_age<<endl;

 //LOCAL_CALCS
 //  cout<<nages<<endl;
 //END_CALCS
   int styr_rec; 
   vector cv_srv1(1,nobs_srv1);      //shelf survey CV
   vector cv_srv2(1,nobs_srv2);      //slope survey CV
   vector cv_srv3(1,nobs_srv3);      //Aleutian Islands survey CV 
 //year
  int i
//age
  int j
//sex
  int k
//
  int ii
  int m
 

 LOCAL_CALCS
   styr_rec=styr-nages+1;
   if(nselages>nages) nselages=nages;
   if(nselages_srv1>nages) nselages_srv1=nages;  
   if(nselages_srv2>nages) nselages_srv2=nages;
   //calculate cv for surveys
    cv_srv1=elem_div(obs_srv1_sd,obs_srv1);   //shelf survey CV
    cv_srv2=elem_div(obs_srv2_sd,obs_srv2);   //slope survey CV
    cv_srv3=elem_div(obs_srv3_sd,obs_srv3);   //Aleutian Island survey CV
   // cout <<cv_srv3<<endl;exit(1);
   //change weights to tons
   wt=wt*.001;
   //  nages=nages-10
   //  cout<<nobs_srv1<<endl;
   //  cout<<wt<<endl;
 END_CALCS

  //3darray obs_p_srv1_age(1,2,1,nobs_srv1_age,1,nages)  //calculated proportion from age bin input
  //3darray obs_p_srv2_age(1,2,1,nobs_srv2_age,1,nages)  // same for slope survey
  //3darray obs_p_srv1_age_r(1,2,1,nobs_srv1_age,1,nages) //calculated with first two ages chopped off
  //3darray obs_p_srv2_age_r(1,2,1,nobs_srv2_age,1,nages) // same for slope survey

  vector obs_sexr(1,nobs_fish_length)  // prop. females in fishery length data
  //vector obs_sexr_srv1(1,nobs_srv1_age) // prop. females in shelf survey ages
  //vector obs_sexr_srv2(1,nobs_srv2_age) // prop. females in slope survey ages
  vector obs_sexr_srv1_2(1,nobs_srv1_length) // prop. males in shelf survey length data
  vector obs_sexr_srv2_2(1,nobs_srv2_length) // prop. males in slope survey length data
  vector obs_sexr_srv3_2(1,nobs_srv3_length) // prop. males in Aleutian Islands length data
  number obs_mean_sexr    //average proportion of males in shelf survey population estimates
  number obs_SD_sexr      //standard deviation from male prop. in shelf survey population estimates
  vector pred_sexr(styr,endyr)   //proportion of males in num at age matrix to be calculated
  
INITIALIZATION_SECTION
  M M_init
  q1 .8   ;
  //can have different mortality for males and females
  F40 .20
  F35 .21
  F30 .23
  // s1_fish .4
  // s1_shelf .4
  // s1_slope .4
  mean_log_rec 10.
  log_avg_fmort -5.
  //proportion of the estimated BSAI survey biomass that is on the shelf
  q2 srv2_qinit   //proportion of the estimated BSAI survey biomass that is on the slope
  q3 .45   //proportion of the estimated BSAI survey biomass that is in the Aleutian Islands
  fmort_dev 0.00001
  fish_slope_f 1.0
  fish_sel50_f  4.
  fish_slope_m  1.0
  fish_sel50_m  4.
  srv1_slope_f1  .8
  srv1_slope_f2  .8
  srv1_slope_m1  .8
  srv1_sel50_f1  2.
  srv1_sel50_f2  8.
  srv1_sel50_m1  2.
  srv1_sel50_m2  8.
  srv1_slope_m2 .4
  srv2_slope_f  .4
  srv2_sel50_f  4.
  srv2_slope_m  .4
  srv2_sel50_m  4.
  srv3_slope_f  .4
  srv3_sel50_f  4.
  srv3_slope_m  .4
  srv3_sel50_m  4.
  alpha 1.
  beta 0.

PARAMETER_SECTION
 //parameters to be estimated are all ones that begin with init_ and have a positive
 //phase, negative phase means are fixed.
 //phase of 8 is greater than last phase so does q1 in last phase  
  // init_bounded_number q1(.5,2,8)
 //fix q1 to be 1 otherwise it went to lower bound of .5
  init_bounded_number q1(0.1,2.0,4)
  init_bounded_number q2(0.01,1.5,-5)
  init_bounded_number q3(0.05,1.5,2)
  init_number alpha(-4)         // used to estimate temperature effect on shelf survey catchability
  init_number beta(phase_Temp)  // used to estimate temperature effect on shelf survey catchability
 //phase of -1 means M is fixed   
  init_bounded_vector M(1,2,.02,.8,phase_M)
  //vector M(1,2)
  init_number mean_log_rec(1)
  // init_bounded_dev_vector rec_dev(styr_rec,endyr-2,-15,15,2)
  init_bounded_dev_vector rec_dev(styr_rec,endyr,-15,15,2) //JNI
//  init_vector rec_dev_future(styr_fut,endyr_fut,phase_F40);
  init_number log_avg_fmort(2)
  !! cout <<log_avg_fmort<<endl;
  init_bounded_dev_vector fmort_dev(styr,endyr,-3,3,1)
 
//  Selectivity parameters 
  init_bounded_number fish_slope_f(.1,5.0,-phase_logistic_sel)
  init_bounded_number fish_sel50_f(1.,25.,phase_logistic_sel+1)
  init_bounded_number fish_slope_m(.05,5.,-phase_logistic_sel)
  init_bounded_number fish_sel50_m(1.,25.,phase_logistic_sel+1)

  // Ascending
  init_bounded_number srv1_slope_f1(.1,5.0,phase_logistic_sel+1)
  init_bounded_number srv1_sel50_f1(1.,10.,phase_logistic_sel)
  // Descending
  init_bounded_number srv1_slope_f2(.1,5.0,phase_logistic_sel+1)
  init_bounded_number srv1_sel50_f2(1,20.,phase_logistic_sel)

  // Ascending
  init_bounded_number srv1_slope_m1(.1,5.,phase_logistic_sel)
  init_bounded_number srv1_sel50_m1(.2,12.,phase_logistic_sel+1)
  // Descending
  init_bounded_number srv1_slope_m2(.1,5.,phase_logistic_sel)
  init_bounded_number srv1_sel50_m2(1.,20.,phase_logistic_sel+1)

  init_bounded_number srv2_slope_f(.1,5.0,phase_logistic_sel)
  init_bounded_number srv2_sel50_f(1.,10.,phase_logistic_sel+1)
  init_bounded_number srv2_slope_m(.1,5.0,phase_logistic_sel)
  init_bounded_number srv2_sel50_m(1.,12.,phase_logistic_sel+1)

  init_bounded_number srv3_slope_f(.01,5.,phase_logistic_sel)
  init_bounded_number srv3_sel50_f(1.,20.,phase_logistic_sel)
  init_bounded_number srv3_slope_m(.01,5.,phase_logistic_sel)
  init_bounded_number srv3_sel50_m(1.,20.,phase_logistic_sel)

// Parameters for computing SPR rates 
  init_bounded_number F40(0.01,1.,phase_F40)
  init_bounded_number F35(0.01,1.,phase_F40)
  init_bounded_number F30(0.01,1.,phase_F40)

//  matrix log_sel_fish(1,2,1,nages)
//  matrix log_sel_srv1(1,2,1,nages)
  matrix sel(1,2,1,nages)
  matrix sel_srv1(1,2,1,nages)
  matrix sel_srv2(1,2,1,nages)
  matrix sel_srv3(1,2,1,nages)
//  vector avgsel_fish(1,2)
  vector avgsel_srv1(1,2)
  vector avgsel_srv2(1,2)
  matrix popn(1,2,styr,endyr)
  matrix totn_srv1(1,2,styr,endyr)
  matrix totn_srv2(1,2,styr,endyr)
  matrix totn_srv3(1,2,styr,endyr)
  vector temp1(1,nages)
  vector temp2(1,nages)
  vector explbiom(styr,endyr)
  vector pred_srv1(styr,endyr)
  vector pred_srv2(styr,endyr)
  vector pred_srv3(styr,endyr)
  3darray pred_p_fish(1,2,1,nobs_fish_length,1,nlen)
//  3darray pred_p_fish_1(1,2,1,styr,endyr,1,nlen)
  3darray pred_p_srv2_age(1,2,1,nobs_srv2_age,1,nages)
  3darray pred_p_srv3_age(1,2,1,nobs_srv3_age,1,nages)
//  3darray pred_p_srv1_age_1(1,2,styr,endyr,1,nages)

  3darray pred_p_srv1_length(1,2,1,nobs_srv1_length,1,nlen)
  3darray pred_p_srv2_length(1,2,1,nobs_srv2_length,1,nlen)
  3darray pred_p_srv3_length(1,2,1,nobs_srv3_length,1,nlen)

  vector pred_catch(styr,endyr)
  3darray natage(1,2,styr,endyr,1,nages)
  3darray catage(1,2,styr,endyr,1,nages)
  //matrix u(styr,endyr,1,nages)
  3darray Z(1,2,styr,endyr,1,nages)
  3darray F(1,2,styr,endyr,1,nages)
  3darray S(1,2,styr,endyr,1,nages)
  vector fmort(styr,endyr)
  number rbar
  vector surv(1,2)  //survival for each sex
  vector offset(1,6)
  number rec_like
  number catch_like
  number sexr_like
  vector age_like(1,2)
  vector length_like(1,4)
  number fpen    
  number surv1_like
  number surv2_like
  number surv3_like
  objective_function_value obj_fun
  number tmp
  vector pred_sexr(styr,endyr)
  vector preds_sexr(styr,endyr)
 // Stuff for SPR and yield projections
  number sigmar
  number ftmp
  number SB0
  number SBF40
  number SBF35
  number SBF30
  number sprpen
  matrix Nspr(1,4,1,nages)
  3darray nage_future(1,2,styr_fut,endyr_fut,1,nages)
  matrix fspbiom_fut(1,4,styr_fut,endyr_fut)
  3darray F_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray Z_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray S_future(1,2,styr_fut,endyr_fut,1,nages)
  3darray catage_future(1,2,styr_fut,endyr_fut,1,nages)
  number avg_rec_dev_future
  vector avg_F_future(1,4)
  sdreport_vector srv_sel_q_f(1,nages);
  sdreport_vector srv_sel_q_m(1,nages);
  sdreport_vector total_biomass(styr,endyr);
  sdreport_vector fspbio(styr,endyr);
  sdreport_number endbiom
  sdreport_number depletion
  sdreport_matrix catch_future(1,3,styr_fut,endyr_fut) // Note, don't project for F=0 (it 
  sdreport_matrix future_biomass(1,4,styr_fut,endyr_fut)
  vector explbiom_fut(styr_fut,endyr_fut)
  number maxsel_fish
  number maxsel_srv1
  number maxsel_srv2
  number maxsel_srv3
  number mlike
  number qlike
  number flike
  vector qtime(styr,endyr)
 
PRELIMINARY_CALCS_SECTION
//compute sex ratio in  catch
  for(i=1; i<=nobs_fish_length;i++)
    obs_sexr(i) = sum(obs_p_fish(1,i))/sum(obs_p_fish(1,i) + obs_p_fish(2,i)); 

//length obs sex ratio in surveys    proportion of males
  for(i=1; i<=nobs_srv1_length;i++)
    obs_sexr_srv1_2(i) = sum(obs_p_srv1_length(2,i))/ sum(obs_p_srv1_length(1,i) + obs_p_srv1_length(2,i));
  for(i=1; i<=nobs_srv2_length;i++)
    obs_sexr_srv2_2(i) = sum(obs_p_srv2_length(2,i))/sum(obs_p_srv2_length(1,i) + obs_p_srv2_length(2,i)); 
  for(i=1; i<=nobs_srv3_length;i++)
    obs_sexr_srv3_2(i) = sum(obs_p_srv3_length(2,i))/sum(obs_p_srv3_length(1,i) + obs_p_srv3_length(2,i)); 

  obs_mean_sexr = mean(obs_sexr_srv1_2);
  obs_SD_sexr   = std_dev(obs_sexr_srv1_2);

 // cout<< " thru sex ratio "<<endl;
 //Compute offset for multinomial and length bin proportions
 // offset is a constant nplog(p) is added to the likelihood     
 // magnitude depends on nsamples(sample size) and p's_
  //k is sex loop
  offset.initialize();
 // fishery length offset and bin proportions 
  for (i=1; i <= nobs_fish_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_fish(1,i)+obs_p_fish(2,i));
    obs_p_fish(1,i) /= sum(obs_p_fish(1,i));
    obs_p_fish(2,i) /= sum(obs_p_fish(2,i));
    for(k=1; k<=2;k++)
      offset(1) -= nsamples_fish(k,i)*obs_p_fish(k,i) * log(obs_p_fish(k,i)+.0001);
  }
 //shelf survey length offset and bin proportions
  for (i=1; i <= nobs_srv1_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv1_length(1,i)+obs_p_srv1_length(2,i));
    obs_p_srv1_length(1,i) = obs_p_srv1_length(1,i) / sumtot; 
    obs_p_srv1_length(2,i) = obs_p_srv1_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(2) -= nsamples_srv1_length(k,i)*obs_p_srv1_length(k,i) * log(obs_p_srv1_length(k,i)+.0001);
  }

//slope survey length offset and bin proportions
  for (i=1; i <= nobs_srv2_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv2_length(1,i)+obs_p_srv2_length(2,i));
    obs_p_srv2_length(1,i) = obs_p_srv2_length(1,i) / sumtot; 
    obs_p_srv2_length(2,i) = obs_p_srv2_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(3) -= nsamples_srv2_length(k,i)*obs_p_srv2_length(k,i) * log(obs_p_srv2_length(k,i)+.0001);
  }

//Aleutian Islands survey length offset and bin proportions
  for (i=1; i <= nobs_srv3_length; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv3_length(1,i)+obs_p_srv3_length(2,i));
    obs_p_srv3_length(1,i) = obs_p_srv3_length(1,i) / sumtot; 
    obs_p_srv3_length(2,i) = obs_p_srv3_length(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(4) -= nsamples_srv3_length(k,i)*obs_p_srv3_length(k,i) * log(obs_p_srv3_length(k,i)+.0001);
  }

//slope survey age offset 
  for (i=1; i <= nobs_srv2_age; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv2_age(1,i)+obs_p_srv2_age(2,i));
    obs_p_srv2_age(1,i) = obs_p_srv2_age(1,i) / sumtot; 
    obs_p_srv2_age(2,i) = obs_p_srv2_age(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(5) -= nsamples_srv2_age(k,i)*obs_p_srv2_age(k,i) * log(obs_p_srv2_age(k,i)+.0001);
  }
//Aleutian Islands survey age offset 
  for (i=1; i <= nobs_srv3_age; i++)
  {
    double sumtot ;
    sumtot = sum(obs_p_srv3_age(1,i)+obs_p_srv3_age(2,i));
    obs_p_srv3_age(1,i) = obs_p_srv3_age(1,i) / sumtot; 
    obs_p_srv3_age(2,i) = obs_p_srv3_age(2,i) / sumtot; 
    for(k=1; k<=2;k++)
      offset(6) -= nsamples_srv3_age(k,i)*obs_p_srv3_age(k,i) * log(obs_p_srv3_age(k,i)+.0001);
  }


PROCEDURE_SECTION
//this is for bootstraping where qrun is a vector of q's from bootstrap irun is the 
//run number.  sets the q (=q1) for the run.
//   q1=qrun(irun);
   get_selectivity();
   get_mortality();
    surv(1)=mfexp(-1.0* M(1));
    surv(2)=mfexp(-1.0* M(2));
   get_numbers_at_age();
   get_catch_at_age();

  if (active(F40))
    compute_spr_rates();
  if (last_phase())
  {
    Future_projections();
  }
  if (sd_phase() || mceval_phase())
   {
     srv_sel_q_f = q1*sel_srv1(1) + q2*sel_srv2(1) + q3*sel_srv3(1) ;
     srv_sel_q_m = q1*sel_srv1(2) + q2*sel_srv2(2) + q3*sel_srv3(2) ;
    if (mceval_phase())
    {
      evalout << obj_fun << " " ;
      // loop over years and print in one long row.
      for (i=styr;i<=endyr;i++)
        evalout<<  fspbio(i) << " " << natage(1,i)*wt(1) + natage(2,i)*wt(2) <<" " << 2*natage(1,i,1) <<" ";
      // hit carriage return on file
      evalout <<  endl;
    }
  }

    evaluate_the_objective_function();
   
FUNCTION get_selectivity
//     logistic selectivity curves, asymptotic for fishery, slope survey and the Aleutian Islands but domed shape for shelf survey 
  for (j=1;j<=nages;j++)
  { 
    sel(1,j)=1./(1.+mfexp(-1.*fish_slope_f*(double(j)-fish_sel50_f)));
    sel(2,j)=1./(1.+mfexp(-1.*fish_slope_m*(double(j)-fish_sel50_m)));
//  ascending limb of curve for shelf survey
    sel_srv1(1,j)=1./(1.+mfexp(-1.*srv1_slope_f1*(double(j)-srv1_sel50_f1)));
    sel_srv1(2,j)=1./(1.+mfexp(-1.*srv1_slope_m1*(double(j)-srv1_sel50_m1)));
//  decending limb of curve for shelf survey 
    temp1=1./(1.+mfexp(srv1_slope_f2*(double(j)-srv1_sel50_f2)));
    temp2=1./(1.+mfexp(srv1_slope_m2*(double(j)-srv1_sel50_m2)));
    sel_srv1(1,j)=sel_srv1(1,j)*temp1(j);
    sel_srv1(2,j)=sel_srv1(2,j)*temp2(j);
//  slope surveys
    sel_srv2(1,j)=1./(1.+mfexp(-1.*srv2_slope_f*(double(j)-srv2_sel50_f)));
    sel_srv2(2,j)=1./(1.+mfexp(-1.*srv2_slope_m*(double(j)-srv2_sel50_m)));
//  Aleutian Islands surveys
    sel_srv3(1,j) = 1./(1.+mfexp(-1.*srv3_slope_f*(double(j)-srv3_sel50_f)));
    sel_srv3(2,j) = 1./(1.+mfexp(-1.*srv3_slope_m*(double(j)-srv3_sel50_m)));
  } 
  // Jim changed this...
  maxsel_fish=1.0; // max(sel(1));
  // if(maxsel_fish<max(sel(2))) maxsel_fish=max(sel(2));
  // obj_fun += 10.*square(log(sel_srv1(1,4)));
  // obj_fun += 10.*square(log(sel_srv1(2,4)));
  sel_srv1(1) /= sel_srv1(1,4); // normalized by 4th age class 
  sel_srv1(2) /= sel_srv1(2,4); // normalized by 4th age class 
  maxsel_srv2=1.0; // max(sel_srv2(1));
  // if(maxsel_srv2<max(sel_srv2(2))) maxsel_srv2=max(sel_srv2(2)); 
  maxsel_srv3=1.0; // max(sel_srv3(1));
  // if(maxsel_srv3<max(sel_srv3(2))) maxsel_srv3=max(sel_srv3(2)); 

FUNCTION get_mortality
  // cout <<log_avg_fmort<<endl<< endl<< mfexp(log_avg_fmort)<<endl<<fmort_dev<<endl;
  fmort = mfexp( log_avg_fmort+fmort_dev);
  for(k=1;k<=2;k++)
  {
    for (i=styr;i<=endyr;i++)
    {
      F(k,i)=sel(k)*fmort(i);
      Z(k,i)=F(k,i) + M(k);
    }
  }
  S = mfexp(-1.*Z);
 // cout<<"to end of get_mortality"<<endl;

FUNCTION get_numbers_at_age
  int itmp;
  //calc initial population  
  for (j=1;j<nages;j++)
  {
    itmp=styr+1-j;
    natage(1,styr,j)=mfexp(mean_log_rec-(M(1)*double(j-1))+rec_dev(itmp));
    natage(2,styr,j)=mfexp(mean_log_rec-(M(2)*double(j-1))+rec_dev(itmp));
    //cout<<"natage"<<natage(1,styr,j)<<endl;
  }
  itmp=styr+1-nages;
  //last age    
  natage(1,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(1)*(nages-1)))/(1.- surv(1));
  natage(2,styr,nages)=mfexp(mean_log_rec+rec_dev(itmp)-(M(2)*(nages-1)))/(1.- surv(2));

 // Now do for next several years----------------------------------
  for (i=styr+1;i<=endyr;i++)
  {
    //for age 1 recruits in the last year use value read in from data file
    natage(1,i,1)=mfexp(mean_log_rec+rec_dev(i));
    natage(2,i,1)=natage(1,i,1);
  }

 //numbers at age
  for(k=1;k<=2;k++)
  {
    for (i=styr;i< endyr;i++)
    {
      //subvector - avoids writing a j loop  =++ increments the right side 
      //(1,nages-1) to 1+1 to nages-1+1 then does the assignment x(i)(1,n) 
      //takes the ith row of x the columns 1 to n
      //      natage(k,i+1)(2,nages)=++elem_prod(natage(k,i)(1,nages-1),S(k,i)(1,nages-1));
      for(j=1;j<nages;j++)
      {
        natage(k,i+1,j+1)=natage(k,i,j)*S(k,i,j); 
      }
      //accumulates oldest ages
      natage(k,i+1,nages)+=natage(k,i,nages)*S(k,i,nages);
      //popn is exploitable numbers
      popn(k,i)= natage(k,i)*sel(k);
    }
    // cout<<"to popn"<<endl; 
    popn(k,endyr)=natage(k,endyr)*sel(k);
  }
  for (i=styr;i<=endyr;i++)
  {
    pred_sexr(i)=sum(natage(2,i))/(sum((natage(1,i)+natage(2,i))));  //calculation of prop. of males in pred. population 
    for(k=1;k<=2;k++)
    {
      totn_srv1(k,i) = q1*(natage(k,i)*sel_srv1(k)); // not used in further calculations
      totn_srv2(k,i) = q2*(natage(k,i)*sel_srv2(k)); // not used in further calculations        
    }
  }
  //predicted survey values
  fspbio.initialize(); 
  explbiom.initialize();
  total_biomass.initialize();
  pred_srv1.initialize();
  pred_srv2.initialize();
  pred_srv3.initialize();
  qtime=q1;
  for (i=styr;i<=endyr;i++)
  {
    fspbio(i) = natage(1,i)*elem_prod(wt(1),maturity);
    if (active(beta))
    {
      if (i>=1991 && i-1990 <= nobs_srv1)      //JNI catchability calculation for survey years
      {
        // qtime(i)  = q1*mfexp(-alpha);
        qtime(i) *= mfexp(beta*bottom_temps(i-yrs_srv1(1)+1));
      }
    }
    for(k=1;k<=2;k++)
    {
      pred_srv1(i) += qtime(i)*(natage(k,i)*elem_prod(sel_srv1(k),wt(k)));   //shelf survey 
      //pred_srv1(i) += q1*(natage(k,i)*elem_prod(sel_srv1(k),wt(k))); //maxsel_srv1;   //shelf survey, without temperature q modeling
      pred_srv2(i) += q2*(natage(k,i)*elem_prod(sel_srv2(k),wt(k))); // /maxsel_srv2;         //slope survey JNI
      pred_srv3(i) += q3*(natage(k,i)*elem_prod(sel_srv3(k),wt(k))); // /maxsel_srv3;         //Aleutian Islands survey JNI

      //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
      // are set different in the tpl file the program will take to value from the bin file and use that 
      //   pred_srv1(i)=1.0*(natage(i)*elem_prod(sel_srv1,wt));
      explbiom(i) += natage(k,i)*elem_prod(sel(k),wt(k))/maxsel_fish;
      total_biomass(i) += natage(k,i)*wt(k);
    }
  }
  // cout <<q3<<endl<<pred_srv3<<endl;exit(1);
  //don't need to divide by max_sel because totn_srv1 is calculated using selectivities and the
  //max_sel would cancel out.

  // Fitting the survey length compositions
  for(i=1; i<=nobs_srv1_length;i++)
  {
    ii = yrs_srv1_length(i);
    pred_p_srv1_length(1,i) = q1 * elem_prod(sel_srv1(1),natage(1,ii)) * lenage(1);
    pred_p_srv1_length(2,i) = q1 * elem_prod(sel_srv1(2),natage(2,ii)) * lenage(2);
    dvariable sum_tot = sum(pred_p_srv1_length(1,i)+pred_p_srv1_length(2,i));
    pred_p_srv1_length(1,i) /= sum_tot;
    pred_p_srv1_length(2,i) /= sum_tot;
  }
 
  for(i=1; i<=nobs_srv2_length;i++)
  {
    ii = yrs_srv2_length(i);
    pred_p_srv2_length(1,i)=q2*elem_prod(sel_srv2(1),natage(1,ii))*lenage(1);
    pred_p_srv2_length(2,i)=q2*elem_prod(sel_srv2(2),natage(2,ii))*lenage(2);
    dvariable sum_tot = sum(pred_p_srv2_length(1,i)+pred_p_srv2_length(2,i));
    pred_p_srv2_length(1,i) /= sum_tot;
    pred_p_srv2_length(2,i) /= sum_tot;
  }
  for(i=1; i<=nobs_srv3_length;i++)
  {
    ii = yrs_srv3_length(i);
    pred_p_srv3_length(1,i)=q3*elem_prod(sel_srv3(1),natage(1,ii))*lenage(1);
    pred_p_srv3_length(2,i)=q3*elem_prod(sel_srv3(2),natage(2,ii))*lenage(2);
    dvariable sum_tot = sum(pred_p_srv3_length(1,i)+pred_p_srv3_length(2,i));
    pred_p_srv3_length(1,i) /= sum_tot;
    pred_p_srv3_length(2,i) /= sum_tot;
  }

  //Calculation of survey age composition
  // BS slope survey 
  for(i=1; i<=nobs_srv2_age;i++)
  {
    ii = yrs_srv2_age(i);
    pred_p_srv2_age(1,i) = q2 * elem_prod(sel_srv2(1),natage(1,ii));
    pred_p_srv2_age(2,i) = q2 * elem_prod(sel_srv2(2),natage(2,ii));
    dvariable sum_tot = sum(pred_p_srv2_age(1,i)+pred_p_srv2_age(2,i));
    pred_p_srv2_age(1,i) /= sum_tot;
    pred_p_srv2_age(2,i) /= sum_tot;
  }
   //  Aleutian Islands survey
  for(i=1; i<=nobs_srv3_age;i++)
  {
    ii = yrs_srv3_age(i);
    pred_p_srv3_age(1,i) = q3 * elem_prod(sel_srv3(1),natage(1,ii));
    pred_p_srv3_age(2,i) = q3 * elem_prod(sel_srv3(2),natage(2,ii));
    dvariable sum_tot = sum(pred_p_srv3_age(1,i)+pred_p_srv3_age(2,i));
    pred_p_srv3_age(1,i) /= sum_tot;
    pred_p_srv3_age(2,i) /= sum_tot;
  }
  depletion=total_biomass(endyr)/total_biomass(styr);
  endbiom=total_biomass(endyr);

FUNCTION get_catch_at_age
  pred_catch.initialize();
  for(k=1;k<=2;k++)
  {      
    for (i=styr; i<=endyr; i++)
    {
      //--Baranov's equation here-----------------------------------
      for (j = 1 ; j<= nages; j++)
      {
        catage(k,i,j)  = natage(k,i,j)*F(k,i,j)*(1.-S(k,i,j))/Z(k,i,j);
        pred_catch(i) += catage(k,i,j)*wt(k,j);
      }
    }
    for (i=1; i<=nobs_fish_length; i++)
    {
      pred_p_fish(k,i)  = elem_prod(sel(k),natage(k,yrs_fish_length(i)))*lenage(k) ;
      pred_p_fish(k,i) /= sum(pred_p_fish(k,i)) ;
    }
  }

FUNCTION Future_projections
  for(k=1;k<=2;k++)
  {
    nage_future(k,styr_fut)(2,nages)=++elem_prod(natage(k,endyr)(1,nages-1),S(k,endyr)(1,nages-1));
    nage_future(k,styr_fut,nages)+=natage(k,endyr,nages)*S(k,endyr,nages);
  }
  future_biomass.initialize();
  catch_future.initialize();
  for (int l=1;l<=4;l++)
  {
    switch (l)
    {
      case 1:
        ftmp=F40;
        break;
      case 2:
        ftmp=F35;
        break;
      case 3:
        ftmp=F30;
        break;
      case 4:
        ftmp.initialize();
        break;
    }
    // Get future F's
    for(k=1;k<=2;k++)
    {
      for (i=endyr+1;i<=endyr_fut;i++)
      {
        for (j=1;j<=nages;j++)
        {
          F_future(k,i,j) = (sel(k,j)/maxsel_fish)*ftmp;
          Z_future(k,i,j) = F_future(k,i,j)+M(k);
          S_future(k,i,j) = exp(-1.*Z_future(k,i,j));
        }
      }
      // Future Recruitment (and spawners)
      for (i=styr_fut;i<endyr_fut;i++)
      {
        nage_future(k,i,1)  = median_rec;
       // Now graduate for the next year....
        nage_future(k,i+1)(2,nages) = ++elem_prod(nage_future(k,i)(1,nages-1),S_future(k,i)(1,nages-1));
        nage_future(k,i+1,nages)   += nage_future(k,i,nages)*S_future(k,i,nages);
      }
      nage_future(k,endyr_fut,1)  = median_rec;
      // Now get catch at future ages
      for (i=styr_fut; i<=endyr_fut; i++)
      {
        for (j = 1 ; j<= nages; j++)
        {
          catage_future(k,i,j) = nage_future(k,i,j) * F_future(k,i,j) * ( 1.- S_future(k,i,j) ) / Z_future(k,i,j);
         if(k==1)
          {
          fspbiom_fut(l,i) += nage_future(1,i,j)*wt(1,j)*maturity(j);
          }
        }
        if (l!=4) catch_future(l,i)   += catage_future(k,i)*wt(k);
        future_biomass(l,i) += nage_future(k,i)*wt(k);
      }   //end loop over future years
    }   //end loop over sex
    fspbiom_fut(l)=0.;
    for(i=styr_fut;i<=endyr_fut;i++)
      fspbiom_fut(l,i) = elem_prod(nage_future(1,i),wt(1)) * maturity;
  }   //End of loop over F's
  
FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0.initialize();
  SBF40.initialize();
  SBF35.initialize();
  SBF30.initialize();

  // Initialize the recruit (1) for each F  (F40 etc)
  for (i=1;i<=3;i++)
    Nspr(i,1)=1.;

  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*exp(-M(1));
    Nspr(2,j)=Nspr(2,j-1)*exp(-(M(1)+F40*sel(1,j-1)/maxsel_fish));
    Nspr(3,j)=Nspr(3,j-1)*exp(-(M(1)+F35*sel(1,j-1)/maxsel_fish));
    Nspr(4,j)=Nspr(4,j-1)*exp(-(M(1)+F30*sel(1,j-1)/maxsel_fish));
  }
 // Now do plus group
  Nspr(1,nages)=Nspr(1,nages-1)*exp(-1.*M(1))/(1.-exp(-1.*M(1)));
  Nspr(2,nages)=Nspr(2,nages-1)*exp(-1.* (M(1)+F40*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F40*sel(1,nages)/maxsel_fish)));
  Nspr(3,nages)=Nspr(3,nages-1)*exp(-1.* (M(1)+F35*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F35*sel(1,nages)/maxsel_fish)));
  Nspr(4,nages)=Nspr(4,nages-1)*exp(-1.* (M(1)+F30*sel(1,nages-1)/maxsel_fish))/ (1.-exp(-1.*(M(1)+F30*sel(1,nages)/maxsel_fish)));
 //cout<<"plus group"<<endl;
  for (j=1;j<=nages;j++)
  {
   // Kill them off till april (0.25) atf spawn in winter so put in 0.0
   //         Number   ProportMat  Wt    Amount die off prior to spawning (within that year)
    SB0    += Nspr(1,j)*maturity(j)*wt(1,j)*mfexp(-0.0*M(1));
    SBF40  += Nspr(2,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F40*sel(1,j)/maxsel_fish));
    SBF35  += Nspr(3,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F35*sel(1,j)/maxsel_fish));
    SBF30  += Nspr(4,j)*maturity(j)*wt(1,j)*mfexp(-0.0*(M(1)+F30*sel(1,j)/maxsel_fish));
  }
  sprpen    = 200.*square((SBF40/SB0)-0.4);
  sprpen   += 200.*square((SBF35/SB0)-0.35);
  sprpen   += 200.*square((SBF30/SB0)-0.30);

FUNCTION evaluate_the_objective_function
  length_like.initialize();
  age_like.initialize();
  fpen.initialize();
  rec_like.initialize();
  surv1_like.initialize();
  surv2_like.initialize();
  surv3_like.initialize();
  catch_like.initialize();
  sexr_like.initialize();

  if (active(rec_dev))
  {
    length_like.initialize();   //length-like vector has the likelihoods for the 4 components: 1) fishery length 2) shelf survey lengths 3) slope survey lengths 4) Aleutians
    int ii;

    //recruitment likelihood - norm2 is sum of square values   
    rec_like = .1* norm2(rec_dev);

    for(k=1;k<=2;k++)
    {
      for (i=1; i <= nobs_fish_length; i++)
      {
        ii=yrs_fish_length(i);
        //fishery length likelihood fitting
        length_like(1) -= nsamples_fish(k,i)*(1e-5+obs_p_fish(k,i))*log(pred_p_fish(k,i)+1e-5);
      }
    }
    //add the offset to the likelihood   
    length_like(1)-=offset(1);

    //shelf survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv1_length; i++)
        length_like(2)-=nsamples_srv1_length(k,i)*(1e-3+obs_p_srv1_length(k,i))*log(pred_p_srv1_length(k,i)+1e-3);
    length_like(2)-=offset(2);

    //slope survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv2_length; i++)
        length_like(3)-=nsamples_srv2_length(k,i)*(1e-3+obs_p_srv2_length(k,i))*log(pred_p_srv2_length(k,i)+1e-3);
    length_like(3)-=offset(3);

    //Aleutian Island survey length composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv3_length; i++)
        length_like(4)-=nsamples_srv3_length(k,i)*(1e-3+obs_p_srv3_length(k,i))*log(pred_p_srv3_length(k,i)+1e-3);
    length_like(4)-=offset(4);

  //slope survey age composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv2_age; i++)
        age_like(1)-=nsamples_srv2_age(k,i)*(1e-3+obs_p_srv2_age(k,i))*log(pred_p_srv2_age(k,i)+1e-3);
    age_like-=offset(5);

 //Aleutian survey age composition fitting
    for(k=1;k<=2;k++)
      for (i=1; i <=nobs_srv3_age; i++)
        age_like(2)-=nsamples_srv3_age(k,i)*(1e-3+obs_p_srv3_age(k,i))*log(pred_p_srv3_age(k,i)+1e-3);
    age_like-=offset(6);
  }
  // Fit to indices (lognormal) 
  //weight each years estimate by 1/(2*variance) - use cv as an approx to s.d. of log(biomass) 

   surv1_like = norm2(elem_div(log(obs_srv1+.01)-log(pred_srv1(yrs_srv1)+.01),sqrt(2)*cv_srv1));
   surv2_like = norm2(elem_div(log(obs_srv2+.01)-log(pred_srv2(yrs_srv2)+.01),sqrt(2)*cv_srv2));
   surv3_like = norm2(elem_div(log(obs_srv3+.01)-log(pred_srv3(yrs_srv3)+.01),sqrt(2)*cv_srv3));
   double var_tmp; for (i=1;i<=nobs_srv3;i++) { var_tmp = 2.*square(log(obs_srv3(i))*cv_srv3(i)); surv3_like += square(log(obs_srv3(i)+.01)-log(pred_srv3(yrs_srv3(i))+.01))/var_tmp; }
   // surv_like = norm2(log(obs_srv1+.01)-log(pred_srv1(yrs_srv1)+.01));

    catch_like=norm2(log(catch_bio+.000001)-log(pred_catch+.000001));

   // sex ratio likelihood
   // This didn't seem correct--why not use all the data?
   // sexr_like=0.5*norm2((obs_mean_sexr-pred_sexr)/obs_SD_sexr); 
   /*
  for(i=1; i<=nobs_fish_length;i++)
   sexr_like += 0.5*square((obs_sexr(i)-pred_sexr(yrs)/obs_SD_sexr); 
    obs_sexr(i) = sum(obs_p_fish(1,i))/sum(obs_p_fish(1,i) + obs_p_fish(2,i)); 
  for(i=1; i<=nobs_srv1_length;i++)
    obs_sexr_srv1_2(i) = sum(obs_p_srv1_length(2,i))/ sum(obs_p_srv1_length(1,i) + obs_p_srv1_length(2,i));
  for(i=1; i<=nobs_srv2_length;i++)
    obs_sexr_srv2_2(i) = sum(obs_p_srv2_length(2,i))/sum(obs_p_srv2_length(1,i) + obs_p_srv2_length(2,i)); 
  for(i=1; i<=nobs_srv3_length;i++)
    obs_sexr_srv3_2(i) = sum(obs_p_srv3_length(2,i))/sum(obs_p_srv3_length(1,i) + obs_p_srv3_length(2,i)); 
   */
 //selectivity likelihood is penalty on how smooth selectivities are   
 //here are taking the sum of squares of the second differences 
 // if(active(log_selcoffs_fish))
 // {  
 //   sel_like(1)=wt_like(1)*norm2(first_difference(first_difference(log_sel_fish(1))));
 //   sel_like(2)=wt_like(2)*norm2(first_difference(first_difference(log_sel_srv1(1))));
 //   sel_like(3)=wt_like(3)*norm2(first_difference(first_difference(log_sel_fish(2))));
 //   sel_like(4)=wt_like(4)*norm2(first_difference(first_difference(log_sel_srv1(2)))); 
 //  for (j=1;j<nages;j++)
 //  {
 //   if(monot_sel==1)
 //   { 
 //    if (log_sel_fish(1,j)>log_sel_fish(1,j+1))
 //       sel_like(1)+=wt_like(5)*square(log_sel_fish(1,j)-log_sel_fish(1,j+1));
 //     if (log_sel_fish(2,j)>log_sel_fish(2,j+1))
 //       sel_like(3)+=wt_like(6)*square(log_sel_fish(2,j)-log_sel_fish(2,j+1));
 //    }
 //   if(monot_sel_srv1==1)
 //    {
 //     if (log_sel_srv1(1,j)>log_sel_srv1(1,j+1))
 //       sel_like(2)+=wt_like(7)*square(log_sel_srv1(1,j)-log_sel_srv1(1,j+1));
 //     if (log_sel_srv1(2,j)>log_sel_srv1(2,j+1))
 //       sel_like(4)+=wt_like(8)*square(log_sel_srv1(2,j)-log_sel_srv1(2,j+1));
 //    }
 //  } 
 //   f+=1.*sum(sel_like);
 //
 //   f+=1.*square(avgsel_fish(1));
 //   f+=1.*square(avgsel_fish(2));
 //   f+=1.*square(avgsel_srv1(1));
 //   f+=1.*square(avgsel_srv1(2));
 // }
  //cout<<" f after catch = "<<f<<endl;
  // Phases less than 3, penalize High F's
    if (current_phase()<2)
    {
       //F's are low for arrowtooth changed the value to compare from .2 to .001
       //don't know if makes any difference since the penalty is reduced at the end
       fpen=10.*norm2(mfexp(fmort_dev+log_avg_fmort)-.01);
    }
    else
    {
      fpen=.001*norm2(mfexp(fmort_dev+log_avg_fmort)-.01);
    }

    if (active(fmort_dev))
    {
      fpen+=.01*norm2(fmort_dev);
    }

  obj_fun += rec_like;
  obj_fun += 1.*length_like(1);     //emphasis factor = 1 for fishery lengths   
  obj_fun += 1.*length_like(2);     //emphasis factor = 1
  obj_fun += 1.*length_like(3);     //emphasis factor = 1
  obj_fun += 1.*length_like(4);     //emphasis factor = 1
  obj_fun += 1.*age_like(1);           //emphasis factor = 1
  obj_fun += 1.*age_like(2);           //emphasis factor = 1
  obj_fun += 1.*surv1_like;         //emphasis factor = 1
  obj_fun += 1.*surv2_like;         //emphasis factor = 1
  obj_fun += 1.*surv3_like;         //emphasis factor = 1
  obj_fun += 300*catch_like;        // large emphasis to fit observed catch
  obj_fun += fpen;
  obj_fun += sprpen;
  obj_fun += 1.*sexr_like;             // male proportion prior, emphasis factor = 1
  /* flike   = 12.5*square(fish_slope_m-fish_slope_f); 
  flike  += 12.5*square(fish_sel50_m-fish_sel50_f); 
  obj_fun += flike;                  // selectivity prior (fishery data lame...
  */ 

//bayesian part normal proirs for M and q1 
//   mlike=mfexp(norm2(M(1)-mfavg)/(2.*(mfcv*mfavg)^2.0))+mfexp(norm2(M(2)-mmavg)/(2.*(mmcv*mmavg)^2.0));
//   qlike=mfexp(norm2(q1-q1avg)/(2.*(q1cv*q1avg)^2.0));

  /*  
  if (active(fish_slope_m))
  {
    flike  = 12.5*square(fish_slope_m-fish_slope_f); 
    flike += 12.5*square(srv1_slope_m2-srv1_slope_f2); 
    flike += 12.5*square(srv2_slope_m-srv2_slope_f); 
    flike += 12.5*square(srv3_slope_m-srv3_slope_f); 
  }*/
//   f=f*mlike*qlike;

FUNCTION write_r
  ivector yrs(styr,endyr);
  yrs.fill_seqadd(styr,1);
  ivector ages(1,nages);
  ages.fill_seqadd(2,1);
  write_R(obj_fun);
  write_R(yrs);
  write_R(ages);
  write_R(maturity);
  write_R(total_biomass);
  dmatrix N_females=value(natage(1));
  dmatrix N_males  =value(natage(2));
  write_R(N_males);
  write_R(N_females);
  dmatrix C_females=value(catage(1));
  dmatrix C_males  =value(catage(2));
  write_R(C_males);
  write_R(C_females);
  dmatrix F_females=value(F(1));
  dmatrix F_males  =value(F(2));
  write_R(F_males);
  write_R(F_females);
  dvector sel_fish_females=value(sel(1));
  dvector sel_fish_males  =value(sel(2));
  write_R(sel_fish_males);
  write_R(sel_fish_females);
  dvector sel_srv1_females=value(sel_srv1(1));
  dvector sel_srv1_males  =value(sel_srv1(2));
  write_R(sel_srv1_males);
  write_R(sel_srv1_females);
  dvector sel_srv2_females=value(sel_srv2(1));
  dvector sel_srv2_males  =value(sel_srv2(2));
  write_R(sel_srv2_males);
  write_R(sel_srv2_females);
  dvector sel_srv3_females=value(sel_srv3(1));
  dvector sel_srv3_males  =value(sel_srv3(2));
  write_R(sel_srv3_males);
  write_R(sel_srv3_females);
  write_R(yrs_srv1);
  write_R(yrs_srv1_length);
  write_R(yrs_srv2_length);
  write_R(yrs_srv3_length);
  write_R(yrs_srv3_age);
  write_R(obs_srv1);
  write_R(obs_srv1_sd);
  write_R(pred_srv1);
  write_R(yrs_srv2);
  write_R(obs_srv2);
  write_R(obs_srv2_sd);
  write_R(pred_srv2);
  write_R(yrs_srv3);
  write_R(obs_srv3);
  write_R(obs_srv3_sd);
  write_R(pred_srv3);
  dmatrix lenage_females = lenage(1);
  dmatrix lenage_males   = lenage(2);
  write_R(lenage_females);
  write_R(lenage_males);

  dmatrix obs_fish_l_females = value( obs_p_fish(1));
  dmatrix obs_fish_l_males   = value( obs_p_fish(2));
  dmatrix pred_fish_l_females = value(pred_p_fish(1));
  dmatrix pred_fish_l_males   = value(pred_p_fish(2));
  write_R(pred_fish_l_females);
  write_R(pred_fish_l_males);
  write_R( obs_fish_l_females);
  write_R( obs_fish_l_males);

  dmatrix obs_srv1_l_females = value( obs_p_srv1_length(1));
  dmatrix obs_srv1_l_males   = value( obs_p_srv1_length(2));
  dmatrix obs_srv2_l_females = value( obs_p_srv2_length(1));
  dmatrix obs_srv2_l_males   = value( obs_p_srv2_length(2));
  dmatrix obs_srv3_l_females = value( obs_p_srv3_length(1));
  dmatrix obs_srv3_l_males   = value( obs_p_srv3_length(2));
  write_R(obs_srv1_l_females);
  write_R(obs_srv1_l_males);
  write_R(obs_srv2_l_females);
  write_R(obs_srv2_l_males);
  write_R(obs_srv3_l_females);
  write_R(obs_srv3_l_males);

  dmatrix pred_srv1_l_females = value(pred_p_srv1_length(1));
  dmatrix pred_srv1_l_males   = value(pred_p_srv1_length(2));
  dmatrix pred_srv2_l_females = value(pred_p_srv2_length(1));
  dmatrix pred_srv2_l_males   = value(pred_p_srv2_length(2));
  dmatrix pred_srv3_l_females = value(pred_p_srv3_length(1));
  dmatrix pred_srv3_l_males   = value(pred_p_srv3_length(2));
  write_R(pred_srv1_l_females);
  write_R(pred_srv1_l_males);
  write_R(pred_srv2_l_females);
  write_R(pred_srv2_l_males);
  write_R(pred_srv3_l_females);
  write_R(pred_srv3_l_males);
  dmatrix obs_srv2_a_females = value( obs_p_srv2_age(1));
  dmatrix obs_srv2_a_males   = value( obs_p_srv2_age(2));
  dmatrix pred_srv2_a_females= value(pred_p_srv2_age(1));
  dmatrix pred_srv2_a_males  = value(pred_p_srv2_age(2));
  dmatrix obs_srv3_a_females = value( obs_p_srv3_age(1));
  dmatrix obs_srv3_a_males   = value( obs_p_srv3_age(2));
  dmatrix pred_srv3_a_females= value(pred_p_srv3_age(1));
  dmatrix pred_srv3_a_males  = value(pred_p_srv3_age(2));
  write_R( obs_srv3_a_females   );
  write_R( obs_srv3_a_males   );
  write_R( pred_srv3_a_females);
  write_R( pred_srv3_a_males  );
  write_R(catch_bio);
  write_R(pred_catch);
  write_R(fspbio);
  write_R(F40);
  write_R(F35);
  write_R(F30);
  write_R(SBF40);
  write_R(SBF35);
  write_R(SBF30);
  write_R(SB0);
  write_R( surv1_like);
  write_R( surv2_like );
  write_R(surv3_like );
  write_R( length_like);
  write_R( rec_like);
  write_R( catch_like);
  write_R( sexr_like);
  write_R( age_like);
  
  write_R( future_biomass );
  write_R(fspbiom_fut );
  write_R( catch_future );
  write_R( q1 );
  write_R( q2 );
  write_R( q3 );
  double M_females = value(M(1));
  double M_males   = value(M(2));
  write_R( M_females);
  write_R( M_males); 
  write_R( pred_sexr );
  write_R( obs_mean_sexr );
  write_R( obs_SD_sexr );

  write_R(alpha);
  write_R( beta );

REPORT_SECTION
  save_gradients(gradients);
  if (last_phase()) write_r();
  report << "Estimated numbers of female fish " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<natage(1,i)<< endl;
   
  report << "Estimated numbers of male fish " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<natage(2,i)<<endl;

  report << "Estimated catch numbers for females " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<catage(1,i)<< endl;
   
  report << "Estimated catch numbers for males  " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<catage(2,i)<<endl;

  report << "Estimated female F mortality " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<F(1,i)<< endl;
   
  report << "Estimated male F mortality " << endl;
    for (i=styr; i<=endyr;i++)
      report <<"   Year: " << i << "," <<F(2,i)<<endl;


  report << "Estimated fishery selectivity for females " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel(1,j)<< endl;
   
  report << "Estimated fishery selectivity for males " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel(2,j)<<endl;

  report << "Estimated shelf survey selectivity for females " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv1(1,j)<< endl;
   
  report << "Estimated shelf survey selectivity for males " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv1(2,j)<<endl;

  report << "Estimated slope survey selectivity for females " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv2(1,j)<< endl;
   
  report << "Estimated slope survey selectivity for males " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv2(2,j)<<endl;


  report << "Estimated Aleutian Islands survey selectivity for females " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv3(1,j)<< endl;
   
  report << "Estimated Aleutian Islands survey selectivity for males " << endl;
    for (j=1; j<=nages;j++)
      report <<"   Age: " << j << "," <<sel_srv3(2,j)<<endl;


  report << endl << "Survey Biomass " << endl;
  report << "Bering Sea shelf survey"  << endl;
  report << "Year" << "," << "observed biomass " << ","<< "predicted biomass" << endl;
    for (i=1; i<=nobs_srv1;i++)
      report << yrs_srv1(i) << ","<< obs_srv1(i) << "," << pred_srv1(yrs_srv1(i)) << endl;
   
  report << "Bering Sea slope survey"  << endl;
  report << "Year" << "," << "observed biomass " << ","<< "predicted biomass" << endl;
    for (i=1; i<=nobs_srv2;i++)
      report << yrs_srv2(i) << ","<< obs_srv2(i) << "," << pred_srv2(yrs_srv2(i)) << endl;
 
  report << "Aleutian Islands survey"  << endl;
  report << "Year" << "," << "observed biomass " << ","<< "predicted biomass" << endl;
    for (i=1; i<=nobs_srv3;i++)
      report << yrs_srv3(i) << ","<< obs_srv3(i) << "," << pred_srv3(yrs_srv3(i)) << endl;

  report <<" Observed female shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << ","<< obs_p_srv1_length(1,i) << endl;

  report <<" Predicted female shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << "," << pred_p_srv1_length(1,i) << endl;

  report <<" Observed male shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << ","<< obs_p_srv1_length(2,i) << endl;

  report <<" Predicted male shelf survey length composition " << endl;
    for (i=1; i<=nobs_srv1_length; i++)
      report << yrs_srv1_length(i) << "," << pred_p_srv1_length(2,i) << endl;
  
  report <<" Observed female slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << ","<< obs_p_srv2_length(1,i) << endl;

  report <<" Predicted female slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << "," << pred_p_srv2_length(1,i) << endl;

  report <<" Observed male slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << ","<< obs_p_srv2_length(2,i) << endl;
  
  report <<" Predicted male slope survey length composition " << endl;
    for (i=1; i<=nobs_srv2_length; i++)
      report << yrs_srv2_length(i) << "," << pred_p_srv2_length(2,i) << endl;
  
  report <<" Observed female Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << ","<< obs_p_srv3_length(1,i) << endl;

  report <<" Predicted female Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << "," << pred_p_srv3_length(1,i) << endl;

  report <<" Observed male Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << ","<< obs_p_srv3_length(2,i) << endl;
  
  report <<" Predicted male Aleutian Islands survey length composition " << endl;
    for (i=1; i<=nobs_srv3_length; i++)
      report << yrs_srv3_length(i) << "," << pred_p_srv3_length(2,i) << endl;

  report <<" Observed female slope survey age composition " << endl;
    for (i=1; i<=nobs_srv2_age; i++)
      report << yrs_srv2_age(i) << ","<< obs_p_srv2_age(1,i) << endl;

  report <<" Predicted female slope survey age composition " << endl;
    for (i=1; i<=nobs_srv2_age; i++)
      report << yrs_srv2_age(i) << "," << pred_p_srv2_age(1,i) << endl;

  report <<" Observed male slope survey age composition " << endl;
    for (i=1; i<=nobs_srv2_age; i++)
      report << yrs_srv2_age(i) << ","<< obs_p_srv2_age(2,i) << endl;

  report <<" Predicted male slope survey age composition " << endl;
    for (i=1; i<=nobs_srv2_age; i++)
      report << yrs_srv2_age(i) << "," << pred_p_srv2_age(2,i) << endl;

  report <<" Observed female Aleutian survey age composition " << endl;
    for (i=1; i<=nobs_srv3_age; i++)
      report << yrs_srv3_age(i) << ","<< obs_p_srv3_age(1,i) << endl;

  report <<" Predicted female Aleutian survey age composition " << endl;
    for (i=1; i<=nobs_srv3_age; i++)
      report << yrs_srv3_age(i) << "," << pred_p_srv3_age(1,i) << endl;

  report <<" Observed male Aleutian survey age composition " << endl;
    for (i=1; i<=nobs_srv3_age; i++)
      report << yrs_srv3_age(i) << ","<< obs_p_srv3_age(2,i) << endl;

  report <<" Predicted male Aleutian survey age composition " << endl;
    for (i=1; i<=nobs_srv3_age; i++)
      report << yrs_srv3_age(i) << "," << pred_p_srv3_age(2,i) << endl;

  report <<" Observed Catch " << endl;
    for (i=styr;  i<=endyr; i++)
       report << "year " << i <<","<< catch_bio(i) << endl;

  report  << "Predicted Catch "<< endl;
    for (i=styr;  i<=endyr;  i++)
        report  <<"year " << i <<","<< pred_catch(i) << endl;

  report  << "Female spawning biomass , Total biomass" <<  endl;
      for (i=styr;  i<=endyr;  i++)
         report << "year " << i <<","<< fspbio(i) <<","<<natage(1,i)*wt(1) + natage(2,i)*wt(2)<< endl;

  report << endl<<endl;
  report << "F40=  " << F40 << endl;
  report << "F35=  " << F35 << endl;
  report << "F30=  " << F30 << endl;
  report << "spawning biomass per recruit at F40 harvest rate "<< SBF40<< endl;
  report << "spawning biomass per recruit at F35 harvest rate "<< SBF35 << endl;
  report << "spawning biomass per recruit at F30 harvest rate "<< SBF30 << endl;
  report << "spawning biomass per recruit with no fishing " << SB0 << endl;

  report << "Likelihood components" << endl;
  report << "shelf survey like component " << surv1_like << endl;
  report << "slope survey like component " << surv2_like << endl;
  report << "Aleutian Islands survey lilke component " <<surv3_like << endl;
  report << "shelf survey length composition " << length_like(2) << endl;
  report << "slope survey length composition " << length_like(3) << endl;
  report << "Aleutian Islands survey length composition " << length_like(4) << endl;
  report << "fishery length composition likelihood " << length_like(1) << endl;
  report << "recruitment likelihood component est.  " << rec_like << endl;
  report << "catch likelihood component est.  "<< catch_like << endl;
  report << "sex ratio likelihood component " << sexr_like << endl;
  report << "slope survey age composition  " << age_like(1) << endl;
  report << "Aleutian Islands survey age composition  " << age_like(2) << endl;
  report << "sprpen " <<sprpen<<endl;
  report << "remainder " <<
  obj_fun - rec_like
          - 1.*length_like(1)
          - 1.*length_like(2)
          - 1.*length_like(3)
          - 1.*length_like(4)
          - 1.*age_like(1)
          - 1.*age_like(2)
          - 1.*surv1_like
          - 1.*surv2_like
          - 1.*surv3_like
          - 300*catch_like
          - fpen
          - sprpen
          - 1.*sexr_like <<endl;             // male proportion prior, emphasis factor = 1

  report << "Projected biomass" << endl;
  report << future_biomass << endl;
  report <<"projected future female spawning biomass " << endl;
  report <<fspbiom_fut << endl;
  report << "Future yield " << endl;
  report << catch_future << endl;
  report << "shelf survey q =" << endl;
  report << q1 << endl;
  report << "slope survey q = " << endl;
  report << q2 << endl;
  report << "Aleutian Islands survey q = " << endl;
  report << q3 << endl;
  report << " female natural mortality for this run" << endl;
  report << M(1) << endl;
  report << " male natural mortality for this run" << endl;
  report << M(2) << endl;

  report <<endl << "temperature effect (q) for the shelf survey "<< endl;
   for (i=1;i<=nobs_srv1;i++)
     report <<yrs_srv1(i)<<","<<bottom_temps(i)<<","<<qtime(yrs_srv1(i))<<endl;

  report << "predicted male proportion in population" << endl;
  for (i=styr;i<=endyr;i++)
    report << i << " " << pred_sexr(i) << endl;

  report << "mean observed prop. male in shelf surveys = "<< obs_mean_sexr << endl;
  report << "stdev of mean observed prop. male in shelf surveys = " << obs_SD_sexr << endl;

  report <<"alpha = "<<alpha<<endl;
  report << "beta= " << beta << endl;
  report<<srv_sel_q_f <<endl;
  report<<srv_sel_q_m <<endl;
  report << " Go drink coffee " << endl;
  //--------------------------------------------------------------------------------------------
  report << "SARA file for Angie Greig" << endl;

  report << "Kamchatka flounder      # stock  " << endl;
  report << "BSAI       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
  report << "2014       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
  report << "1a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
  report << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
  report << "full       # UPDATE (new benchmark full partial)" << endl;
  report << "2          # LIFE_HIST - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "3          # ASSES_FREQ - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "4          # ASSES_LEV - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "5          # CATCH_DAT - SAIP ratings (0 1 2 3 4 5) " << endl;
  report << "3          # ABUND_DAT - SAIP ratings (0 1 2 3 4 5)" << endl;
  report << "inputinNov # Minimum B  Lower 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "input      # Maximum B  Upper 95% confidence interval for spawning biomass in assessment year" << endl;
  report << "input      # BMSY  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
  report << "ADMB       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
  report << "NA         # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
  report << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
  report << "1          # number of fisheries" << endl;
  report << "1000       # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
  report << "1          # recruitment age used by model or size" << endl;
  report << "1          # age+ or mmCW+ used for biomass estimate" << endl;
  report << "\"Single age\"        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
  report << "\"Age model\"         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
  report << "\"Age of maximum F\"  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
  report << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
  report << "ALL" << endl; 

  report <<"#FISHERYYEAR - list years used in the model " << endl;
   for (i=styr;  i<=endyr; i++)
      report << i << "	";
      report<<endl;  

  report<<"#AGE - list of ages used in the model"<<endl;
   for (i=2; i<=25;i++)
      report << i << "	";
      report<<endl; 

  report <<"#RECRUITMENT - Number of recruits by year " << endl;
   for (i=styr;  i<=endyr;  i++)
	   report  << natage(1,i,1)+natage(2,i,1) << "	";
	   report<<endl;     
	
  report <<"#SPAWNBIOMASS - Spawning biomass by year in metric tons " << endl;
   for (i=styr;  i<=endyr;  i++)
      report  << fspbio(i) << "	";
      report<<endl;  

  report <<"#TOTALBIOMASS - Total biomass by year in metric tons " << endl;
   for (i=styr;  i<=endyr;  i++)
      report  << natage(1,i)*wt(1)+natage(2,i)*wt(2) << "	";
      report<<endl;
	
  report <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
	for (i=styr;  i<=endyr;  i++)
	   report  << (F(1,i,24)+ F(2,i,24))/2<< "	";
	   report<<endl;
	  
  report <<"#TOTALCATCH - Total catch by year in metric tons " << endl;
   for (i=styr;  i<=endyr;  i++)
      report  << catch_bio(i) << "	";
      report<<endl;
  
  report <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
      for (i=1;  i<=24;  i++) 
      report  << maturity(i) <<"	";
      report<< endl; 

  report <<"#SPAWNWT - Average spawning weight (in kg) for ages 7-24"<< endl; 
      report <<"0.614	0.811	1.0244	1.251	1.487	1.728	1.971	2.213	2.451	2.684	2.910	3.128	3.337	3.535	3.724	3.902	4.070	4.227	4.375";
      report<<endl;

  report <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
     for (i=1;  i<=24;  i++) 
  report  << 0.11 <<"	";
  report<< endl;   
    for (i=1;  i<=24;  i++) 
  report  << 0.11 <<"	";
  report<< endl;

  report << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
    for (i=styr; i<=endyr;i++)
    report <<natage(1,i)<< "	";
    report<<endl;

  for (i=styr; i<=endyr;i++)
    report <<natage(2,i)<< "	";
    report<<endl;

  report <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
    for (i=1;  i<=24; i++)   
  report <<wt(1,i) << "	";
  report<<endl;
         
    for (i=1; i<=24; i++)
  report<<wt(2,i) << " ";
  report<<endl; 

  report << "#SELECTIVITY - Estimated fishery selectivity for females (first) males (second) at age " << endl;
    for (j=1;  j<=nages;  j++)
     report <<" "<<sel(1,j) << "       ";
     report<<endl;  

    for (j=1;  j<=nages;  j++)
     report <<" "<<sel(2,j)<< "      ";
     report<<endl;
  
  report << "#SURVEYDESC"<<endl;
  report<<"EBS_trawl_survey "<<endl;

  report<<"#SURVEYMULT"<<endl;
  report<<"1"<<endl;

  report << "#EBS_trawl_survey - Bering Sea shelf survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
   for (i=1; i<=nobs_srv1; i++)
     report << yrs_srv1(i) << "	";
     report<<endl;
   for (i=1; i<=nobs_srv1; i++) 
     report<< obs_srv1(i)<< "	";
     report<< endl;
  
  report << "#BS_slope_trawl_survey - Bering Sea slope survey biomass (Year, Obs_biomass, Pred_biomass) " << endl;
   for (i=1; i<=nobs_srv2;i++)
     report << yrs_srv2(i) << "	";
     report<<endl;
   for (i=1; i<=nobs_srv2;i++)
     report << obs_srv2(i) << "	";
     report<<endl;

  report << "#AI_trawl_survey - Aleutian Islands survey biomass (Year, Obs_biomass, Pred_biomass) "  << endl;
   for (i=1; i<=nobs_srv3;i++)
     report << yrs_srv3(i) << "	"; 
     report<<endl;
   for (i=1; i<=nobs_srv3;i++)
     report << obs_srv3(i) << "	";
     report<<endl;	 

  report<<"#STOCKNOTES"<<endl;
  report<<"\"SAFE report indicates that this stock was not subjected to overfishing in 2013 and is neither overfished nor approaching a condition of being overfished in 2014.\""<<endl;

GLOBALS_SECTION
  #include <admodel.h>
	#undef REPORT
	#define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("Kam_R.rep");

RUNTIME_SECTION
  maximum_function_evaluations 4000
  convergence_criteria 1e-3 1e-4 1e-7

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(300);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);

FINAL_SECTION
 // write out sderport stuff
  dvector fspbio_lb=value(elem_div(fspbio,exp(2.*sqrt(log(1+elem_div(square(fspbio.sd),square(fspbio)))))));
  dvector fspbio_ub=value(elem_prod(fspbio,exp(2.*sqrt(log(1+elem_div(square(fspbio.sd),square(fspbio)))))));
  write_R(fspbio_lb);
  write_R(fspbio_ub);
  dvector total_biomass_lb=value(elem_div(total_biomass,exp(2.*sqrt(log(1+elem_div(square(total_biomass.sd),square(total_biomass)))))));
  dvector total_biomass_ub=value(elem_prod(total_biomass,exp(2.*sqrt(log(1+elem_div(square(total_biomass.sd),square(total_biomass)))))));
  write_R(total_biomass_lb);
  write_R(total_biomass);
  write_R(total_biomass_ub);
  write_R(total_biomass.sd);
