Steps and notes for running GOA dusky rockfish assessment
Updated: November 2024

2022 full assessment completed by Ben Williams.
2024 full assessment completed by Kristen Omori
Most steps are completed in the R script yyyy_analysis.R using data queries and
organization from afscassess, but currently ADMB model runs with changes run separately

1. Create new assessment year folder, copy over code folder from previous full assessment

2. In R script, yyyy_analysis.R, update:
   A. Global values (year, TAC, previously accepted model)

3. Run R script to setup data and model.
   A. Verify globals are correct (rec_age, plus_age, hold off on len_bins)
   B. afscassess::sp_switch(species): verify using correct species
   C. afscdata::setup_folders(year): set up data folder structure (raw, user_input, output, sara, sql)
      Run the setup_folders only once
   D. Transfer necessary user_input files into user_input file folder
      i. copy paste all previous years necessary user_input files
	  (len_bins, abc_tac (manually updated), goa_dusk_catch_1977_1990.csv, reader_tester, vast model...)
	  ii. update necessary files (abc_tac, reader_tester, updated vast model)
   E. afscdata::accepted_model(base_year= year-2, base_model= prev_model_name, year= year)
      This will copy over the previously accepted model from the last full assessment for comparison
   F. afscassess::setup_tpl(year): Transfer over associated ADMB file (do some by hand)
      i. data > models > ... make sure the .tpl transfers
	  ii. copy over .exe

4. Pull data (data queries) and data organization/ restructure for model.
   A. goa_dusk(year, off_yr= FALSE): will run all the raw data queries for dusky data > raw
   B. Format all raw data queries. Processed data located in: data > output
   C. Note: survey biomass index now using VAST, but good to pull GAP survey biomass + rema for comparison
   D. Do data comparison from previous year to new data year to ensure no major data changes and queried correctly
   E. Setup .ctl file with new data

5. Setup initial assessment model run
   A. In new labeled model folder (year > model_name_folder, e.g., dusk_test_2>2022), copy over base model
   B. Make sure the .exe, .tpl, .ctl, .data (in .cpp), .pin (for starting values) are copied into model folder
      (delete other output files)
   C. In .ctl file: change the following lines
      i. Line 2. Data file name
	  ii. Line 5. endyr
	  iii. Line 59. yield_ratio (update to output 'yld' yield_ratio
   D. In .tpl file: change the following line
      i. Line 24. change the file name to new file name from the .tpl ## Changing this .ctl = goa_dusk.ctl
	     Now don't need to compile model again
   E. Run model to make sure it runs with outputs in ADMB
   F. Run model directly in ADMB (in shell):
      i. Open ADMB in the 'Command Prompt' [or run_admb('base', verbose=TRUE) or shell('base.exe') in R]
	      a. if in ADMB:
		     1) cd C:\Users\Kristen.Omori\Work\Assessments\...full path...\m22.3a --> this will change directory/ path
             2) type: admb base.tpl  ---> this will make the .exe file
             3) type: base.exe  ---> this will run the admb model (make sure it runs!)
	  ii. Prompt lines:
		  a. admb modelname  (i.e., base.exe)
   G. Can take initial run par file and use for .pin file to see if results change (run without MCMC) & to help convergence
         * copy .par file and rename with .pin (acts as starting values for parameters)
   H. Check if gradients all converged, parameters aren't against bounds.
   I. Resave new .par file as new .pin file (keep old pin files) and rerun model to ensure best model fit until stablized
6. Run process_results() "Process Results" without MCMC/ projections to examine results
   A. process_results() --> reads in output files and format into readable files
   B. output in: year>model>processed
   C. Plot data fits with estimates (survey estimates with data input, biomass comparison, age/length comp estimates with input data)
   D. Once you have good convergence this becomes the final 'base model' for this year's SAFE.

7. Run MCMC
   A. Create new mcmc inside the model folder, and copy over .exe, .ctl, .tpl, and .dat files.
   B. Run mcmc (see above for running prompt lines in ADMB)
      i. with -mcmc 10000000 -mcsave 2000
	  ii. Prompt lines example:
		  a. admb modelname
		  b. model name -mcmc 100000 -mcsave 2000
		  c. model name -mceval
	  iii. the -mceval args tells ADMB which parameters to store/ track for the MCMC
	  iv. output from -mceval = 'evalout.prj' (can rename as evalout.dat)
   C. from D: Once MCMCs are run, open tem.tpl and run with args -mceval (there is a line in the tpl telling
	              ADMB which parameters to store/track for MCMC).
	 D. May not converge, use pin files and/or fix sigma_r

8. Run process_results() function WITH mcmc comparisons or projections to examine results
   A. pulls in results and prcoesses them; output in: year>model>processed
9. Plot the rest of the plots to ensure fit with mcmc:
   A. afscassess::base_plots() can plot necessary plots, can also do individaully

10. Run Retrospective Analysis
   A. run_retro_dusk(year, model = "m22.5a", tpl_name = "base", n_retro = 10, mcmc = 50000000, mcsave = 2000)
   B. process_retro_dusk(). Process retrospective results
   C. Comparison plots for retrospective analysis

11. Projections:
   A. afscassess::run_proj()  ** or can use Jim's projection model

12. Run best_f() function for SARA files:
   A. see Pete's pop... 257 https://github.com/pete-hulson/goa_pop/blob/main/2023/R/2023_analysis.R

13. Plot comparison, model comparison, and model sensitivity/ diagnostics:
   A. Look at likelihood profiles
   B. Dan looks at step-wise data sensitivity by adding data sets in one by one and compares models

14. Run apportionment
   A. Requires download of egoa fractions from AKFIN (format as CSV w/o header material; colnames unchanged)
   B. Using VAST or db-index + rema for Western, Central, Eastern biomass proportions
   C. Applying EGOA fractions to EGOA for further split

15. Final tables and figures:
   A. Run all final tables (check formatting in Rmd)
   B. Run all figures; check folder links in Rmd
   C. Other tables: fishery catch outside of area, bycatch...
   D. Appendix other non-commercial catch

16. SAFE writing and Rmarkdown document.
   A. Make sure to update goa_dusky.Rmd (and can do html for quicker viewing)
   B. Final doc, use word. Check numbers and values match and figure/ table labels are correct.

17. Other associated assessment tasks:
   A. SARA file verify numbers
   B. Director's briefing presentation (2 slides in .ppt). Use previous year's template
   C. Update github with final models, etc...

18. Post assessment
   A. Survey ageing requests in ~ Feb. in AGP portal
   B. Fishery observer special requests ~ Aug. via email
