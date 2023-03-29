###################################################################################################
#  Updates of the RTQC and adjustment procedures for the CHLA parameter 
#
#
# Inputs : FLUORESCENCE_CHLA, CHLA from a B file 
#  
# Outputs : 
#           PARAMETERS : 
#		CHLA_FLUORESCENCE
#		CHLA_FLUORESCENCE_ADJUSTED
#		CHLA_ADJUSTED
# 	    QC : 
#		CHLA_QC 
#		CHLA_FLUORESCENCE_QC 
#		CHLA_ADJUSTED_QC 
#		CHLA_FLUORESCENCE_ADJUSTED_QC
#
#
########################################################################################

library(stringr)

library(ncdf4)

source("./MLD_calc.R") # estimation of the MLD with CTD profile

source("./RunningFilter.R") # median filter
 
source("./julian_to_asol.R") # estimation of the zenith angle according to the LATITUDE, LONGITUDE, JULD of the PROFILE

source("./read_META_CHLA.R") # read SCALE_CHLA and FACTORY_DARK_CHLA in meta file 

source("./DARK_CHLA.R") # DARK work 

source("./NPQ_CHLA.R") # Quenching correction 

uf=commandArgs()

mission  <- uf[2]

##################################################################
##  GLOBAL RANGE 
##################################################################

MIN_RANGE_CHLA=-0.2

MAX_RANGE_CHLA=100.

####################################################################
### Read the PREDEPLOYMENT_CALIB_COEFFICIENT in the metadata file
####################################################################

meta_filename=paste("../DATA/",mission,"/",mission,"_meta.nc",sep="")

filenc_meta=nc_open(meta_filename,readunlim=FALSE,write=FALSE)

CALIB_CHLA=read_META_CHLA(filenc_meta)

FACTORY_DARK_CHLA=CALIB_CHLA$DARK

FACTORY_SCALE_CHLA=CALIB_CHLA$SCALE

####################################################################
## Initialize i_count_dark (it should reach 5 ) 
####################################################################

icount_dark=0 

PRELIM_DARK=rep(NA,5)

DARK_PRIM_CHLA=FACTORY_DARK_CHLA ### the DARK_PRIM_CHLA to apply is initialize at the factory calibration dark 

DARK_PRIM_CHLA_QC="2"  ### even with no npq correction" 

####################################################################
#### Creating the list of files for which we need to recompute  
####################################################################

liste_to_do=read.table("./liste_all_B",header=FALSE, as.is=TRUE)

# List of the file to process
LIST_nc=liste_to_do$V1

### loop on all files of the mission 
for (IDnc in LIST_nc) {

### IDnc 

	print(IDnc)

#### Getting the name of the core file
	file_in_C=str_replace(IDnc,"/B","/")

# if B and C are not in the same mode 
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/R","profiles/D")
	if (!file.exists(file_in_C)) file_in_C=str_replace(file_in_C,"profiles/D","profiles/R")

#########################
# Open the nc File
#########################

# Open the C file
	filenc_C=nc_open(file_in_C,readunlim=FALSE,write=FALSE)

# Open the B file
	filenc_B=nc_open(IDnc,readunlim=FALSE,write=TRUE)


############################################################
## work on index 
#############################################################

#### Get the list of parameters in the profile

	STATION_PARAMETERS=ncvar_get(filenc_B,"STATION_PARAMETERS")

# Stations parameters has a fixed length 64 characters 

	CHLA_STRING=str_pad("CHLA",64,"right")

# Find the profile containing CHLA/BBP and PAR  

	index_chla=which(STATION_PARAMETERS == CHLA_STRING, arr.ind=TRUE)


	if (length(index_chla)==0) {

		next # jump to the next profile if no chla in the profile

	} else {

		iprof_chla=index_chla[,2]

	} 

##########################################################
# MLD estimation
##########################################################

#### Read the C file to estimate the MLD

	TEMP_CTD=ncvar_get(filenc_C,"TEMP")

	PSAL_CTD=ncvar_get(filenc_C,"PSAL")

	PRES_CTD=ncvar_get(filenc_C,"PRES")

	MLD=MLD_calc(PRES_CTD, PSAL_CTD , TEMP_CTD)

	print(MLD)


##########################################################
### Estimate the size of the sliding median 
##########################################################
	
 	DELTA_PRES=diff(PRES_CTD[,iprof_chla])

	PROF_RES=median(DELTA_PRES,na.rm=TRUE)

	if(PROF_RES <=1.) median_window=5

	if(PROF_RES>1 & PROF_RES<3.) median_window=3

 	if(PROF_RES>=3) median_window=2

#############################################################################
# SOLAR angle estimation to light on/off the quenching correction
#############################################################################

	JULD=unique(ncvar_get(filenc_B,"JULD"))
	
	LATITUDE=unique(ncvar_get(filenc_B,"LATITUDE"))

	LONGITUDE=unique(ncvar_get(filenc_B,"LONGITUDE"))

	ASOL=solar_angle_test(JULD,LATITUDE,LONGITUDE)

	FLAG_QUENCHING=ASOL$BOOL

	ASOL=ASOL$ASOL

#	print(FLAG_QUENCHING)

#############################################################################
# READ the CHLA variable and estimate the number of vertical levels
#############################################################################

	CHLA=ncvar_get(filenc_B,"CHLA")

	CHLA_QC=ncvar_get(filenc_B,"CHLA_QC")

	CHLA_ADJUSTED=ncvar_get(filenc_B,"CHLA_ADJUSTED")

	CHLA_ADJUSTED_QC=ncvar_get(filenc_B,"CHLA_ADJUSTED_QC")

	FLUORESCENCE_CHLA=ncvar_get(filenc_B,"FLUORESCENCE_CHLA")

	### PRESSURE FOR CHLA

	PRES_CHLA=PRES_CTD[,iprof_chla]

	FLUORESCENCE=FLUORESCENCE_CHLA[,iprof_chla]
	
	#### Number of measurements 

	Ndepth=length(PRES_CHLA[!is.na(PRES_CHLA)])

###########################################################
##      DARK WORK 
###########################################################

	if ( max(PRES_CHLA,na.rm=TRUE) >= 950 & icount_dark < 5 ) {

		icount_dark=icount_dark+1  ### Test 3 ###

		iDARK=dark_chla(FLUORESCENCE,median_window)

		PRELIM_DARK[icount_dark] = iDARK # Test4 

#		print(PRELIM_DARK)

		DARK_PRIM_CHLA=median(PRELIM_DARK,na.rm=TRUE)

		if ( icount_dark == 5 ) {

			if ( abs( DARK_PRIM_CHLA - FACTORY_DARK_CHLA) >= 0.25 * FACTORY_DARK_CHLA) {

				DARK_PRIM_CHLA_QC="3"   ### if float_dark_chla is too far from calibration there is an issue 

			} else {

				DARK_PRIM_CHLA_QC="1"			

			}

		}
	}

############################################################
##  Apply the DARK Value to the adjusted Field
##############################################################

	CHLA_ADJUSTED=FACTORY_SCALE_CHLA*(FLUORESCENCE_CHLA-DARK_PRIM_CHLA)/2

######  AND DEFINE CHLA_FLUORESCENCE PARAMETER 

	CHLA_FLUORESCENCE = CHLA

	CHLA_FLUORESCENCE_ADJUSTED = 2*CHLA_ADJUSTED #### no reason to multiply by 2 , this is fluorescence !!!! be careful of the slope

###########################################################
##  Init the QC  
###########################################################

	CHLA_QC_value = rep("3", length(PRES_CHLA) )

	if (FLAG_QUENCHING & MLD=0) {   # sun and no way to determine MLD MLD =0  not very happy with it 

		CHLA_ADJUSTED_QC_value = rep("3", length(PRES_CHLA) )   

	} else { 
	
		CHLA_ADJUSTED_QC_value = rep(DARK_PRIM_CHLA_QC, length(PRES_CHLA) )

	}

	CHLA_FLUORESCENCE_QC_value = rep("1", length(PRES_CHLA) )  

	CHLA_FLUORESCENCE_ADJUSTED_QC_value = rep("1", length(PRES_CHLA) )

#############################################################
##	# sun and no way to determine MLD 
#############################################################

	if (FLAG_QUENCHING & MLD=0) CHLA_ADJUSTED_QC_value = rep("3", length(PRES_CHLA) )   

###########################################################
# GLOBAL RANGE TEST 
###########################################################

	for(p in 1:Ndepth){

		if ( (CHLA[p,iprof_chla] < MIN_RANGE_CHLA) | (CHLA[p,iprof_chla] > MAX_RANGE_CHLA) )  {

			CHLA_QC_value[p]=4

			CHLA_ADJUSTED_QC_value[p]=4

			CHLA_FLUORESCENCE_ADJUSTED_QC_value[p]=4

			CHLA_FLUORESCENCE_QC_value[p]=4

		}

	}	

#############################################################
# Quenching Correction 
#############################################################
	if (FLAG_QUENCHING & MLD>0) {

		NPQ=npq_chla(PRES_CHLA,CHLA_ADJUSTED[,iprof_chla],MLD,median_window)

		for (p in 1:NPQ$INDEX){

			CHLA_ADJUSTED[p,iprof_chla]=NPQ$VALUE

			CHLA_ADJUSTED_QC_value[p]=5

		}

	}


#########################################
# PLOT
#########################################

	MINDEPTH=max(PRES_CHLA)
	MAXCHLA=2*max(CHLA_ADJUSTED,na.rm=TRUE)
	#MINDEPTH=MLD+50
	MINDEPTH=400

	name_file=str_sub(IDnc,str_length(IDnc)-15,str_length(IDnc)-3)

	png(file=paste(name_file,".png",sep=""))

	matplot(2*CHLA_ADJUSTED[,iprof_chla],PRES_CHLA,col=2,lwd=2,type="l",ylab="Depth [m]",cex.lab=1.5,cex.axis=1.5,xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, MINDEPTH)))
	matplot(CHLA[,iprof_chla],PRES_CHLA,col=5,lwd=2,type="l",ylab="Depth [m]",cex.lab=1.5,cex.axis=1.5,xlab=expression("Chlorophyll a [mg."*m ^ -3 * "]"),xlim=c(-0.2,MAXCHLA+0.5),ylim=rev(c(0, MINDEPTH)),add=TRUE)

	legend("bottomright",c("CHLA_NPQ","CHLA"),pch=c(".","."),lwd=c(2,2),col=c(2,5),lty=c(1,1),cex=1.2)


	dev.off()

###########################################################################
#	Create the new variables in the nc file
###########################################################################



############################################################################
### Write the DATA in the file 
############################################################################
####  DATA in the file

####  QC in the file 

	CHLA_QC[iprof_chla]=paste(CHLA_QC_value,collapse="")

	CHLA_ADJUSTED_QC[iprof_chla]=paste(CHLA_ADJUSTED_QC_value,collapse="")

#	print(CHLA_QC)

#	print(CHLA_ADJUSTED_QC[iprof_chla])

	ncvar_put(filenc_B,"CHLA_QC",CHLA_QC)

	ncvar_put(filenc_B,"CHLA_ADJUSTED_QC",CHLA_ADJUSTED_QC)

	ncvar_put(filenc_B,"CHLA_ADJUSTED",CHLA_ADJUSTED)

	QC="Done and Write"


	nc_close(filenc_C) # closing C file
	
	nc_close(filenc_B) # closing B file 


} # end loop on all files 

nc_close(filenc_meta) # closing metadata file
