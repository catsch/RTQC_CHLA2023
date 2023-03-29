dark_chla <- function(FLUORESCENCE_CHLA,median_window) {

########### Filter the data 

MED_FLUO_CHLA=RunningFilter(median_window,FLUORESCENCE_CHLA,na.fill=T, ends.fill=T, Method="Median") # From Equation 5 

#### ALT 1 # RES=FLUORESCENCE_CHLA-MED_FLUO_CHLA ## Equation 4

#### ALT 1 # Q10=rep(quantile(RES,0.10,na.rm=TRUE),length(FLUORESCENCE_CHLA))

############ Remove from the estimation of the dark outliers

#### ALT 1 # FLUORESCENCE_CHLA[RES<Q10]=NA #### Remove outliers (not sure that is the best way )

############ Estimate the minimum of the profile without outliers

iDARK=min(FLUORESCENCE_CHLA,na.rm=TRUE) ## Equation 6

return(iDARK)

}
