npq_chla <- function(PRES,CHLA_ADJUSTED,MLD,median_window) {

########### Filter the data 

MED_CHLA_ADJUSTED=RunningFilter(median_window,CHLA_ADJUSTED,na.fill=T, ends.fill=T, Method="Median") # From Equation 5 

RES=CHLA_ADJUSTED-MED_CHLA_ADJUSTED ## Equation 4

Q90=rep(quantile(RES,0.90,na.rm=TRUE),length(CHLA_ADJUSTED))

############ Remove from the estimation of the dark outliers

CHLA_ADJUSTED[RES>Q90]=NA #### Remove outliers (not sure that is the best way )

############ Estimate the minimum of the profile without outliers

print(MLD)

#print(CHLA_ADJUSTED[PRES>=MLD])

#print(PRES)

npq_value=max(CHLA_ADJUSTED[PRES<=0.9*MLD],na.rm=TRUE)
#print(npq_value)

npq_index=which.max(CHLA_ADJUSTED[PRES<=0.9*MLD])


return(list("INDEX"= npq_index,"VALUE"=npq_value))

}
