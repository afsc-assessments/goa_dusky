Base model from 2022 with updated lognormal survey biomass error structure and predicted recruitment start year
in the population projection component inside assessment model.

1. In .ctl file (goa_dr_2022.ctl), change 'ph_logsurv' (line 54) from 1 (normal) to 2 (lognormal).
2. In the 'Get_Population_Projection' function, the starting number at age values for projections...
   Base: pred_rec(1979,endyr-recage)
   Updated model: pred_rec(1977+recage, endyr-recage)