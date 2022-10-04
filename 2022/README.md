# Gulf of Alaska Dusky Rockfish  ASsessment notes

## Models 

 - base = m15.5a w/updated 2022 data 
 - m22.1 = increased length plus group - GAP default VAST
 - m22.2 = increased age plus group - GAP default VAST
 - m22.3 = increased age and length plus groups - GAP default VAST
 - m22.3a = increased age and length plus groups - lognormal error VAST
 

## Model run notes:

The GOA dusky model estimates sigma r.
When running mcmc need to get the sigr estimate from the .par file, paste it into the .ctl file and turn the phase to -1 then run the mcmc.
Todo: diagnose the model and figure out why this is occurring, possibly sigr isn't well defined and maybe best not to estimate it? (Is  not estimated in the northern rockfish model)
