npq_chla <- function(PRES,CHLA_ADJUSTED,MLD,median_window) {

########### Filter the data 

MED_CHLA_ADJUSTED=RunningFilter(median_window,CHLA_ADJUSTED,na.fill=T, ends.fill=T, Method="Median") # From Equation 12 

#####  ALT 1 # RES=CHLA_ADJUSTED-MED_CHLA_ADJUSTED ## Equation 4

#####  ALT 1 # Q90=rep(quantile(RES,0.90,na.rm=TRUE),length(CHLA_ADJUSTED))

############ Remove from the estimation of the dark outliers

#####  ALT 1 # CHLA_ADJUSTED[RES>Q90]=NA #### Remove outliers (not sure that is the best way )

############ Estimate the minimum of the profile without outliers

print(MLD)

#print(CHLA_ADJUSTED[PRES>=MLD])

#print(PRES)

#####  ALT 1 # npq_value=max(CHLA_ADJUSTED[PRES<=0.9*MLD],na.rm=TRUE)

#####  ALT 2 #

npq_value=max(MED_CHLA_ADJUSTED[PRES<=0.9*MLD],na.rm=TRUE)   # Equation 14

#print(npq_value)

####   ALT 1 # npq_index=which.max(CHLA_ADJUSTED[PRES<=0.9*MLD])

#####  ALT 2 #

npq_index=which.max(MED_CHLA_ADJUSTED[PRES<=0.9*MLD])

return(list("INDEX"= npq_index,"VALUE"=npq_value))

}
