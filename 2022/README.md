# Gulf of Alaska Dusky Rockfish Assessment notes

## Models 
 
 - db = design-based trawl survey (1990-2021), data through 2022
 - db.1 = db w/increased length plus group
 - db.2 = db w/increased age plus group
 - db.3 = db w/increased age & length plus groups
 - m15.5a = base w/updated 2022 data 
 - m22 = GAP default VAST (1990+)
 - m22a = lognormal error VAST (1990+)
 - m22.1 = increased length plus group - GAP default VAST
 - m22.1a = increased length plus group - lognormal error VAST
 - m22.2 = increased age plus group - GAP default VAST
 - m22.2a = increased age plus group - lognormal error VAST
 - m22.3 = increased age and length plus groups - GAP default VAST
 - m22.3a = increased age and length plus groups - lognormal error VAST
 
 *Note: Had to fix maturity params for m15.5a otherwise they changed to unrealistic values*

## Model run notes:

The GOA dusky model estimates sigma r.
When running mcmc need to get the sigr estimate from the .par file, paste it into the .ctl file and turn the phase to -1 then run the mcmc.
Todo: diagnose the model and figure out why this is occurring, possibly sigr isn't well defined and maybe best not to estimate it? (Is not estimated in the northern rockfish model)

## Future

Be certain GAP removes 1980s survey data.  
Explore reweighting comp data  
Incorporate survey input sample size from swo work  


