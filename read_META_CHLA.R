read_META_CHLA <- function(filenc_meta) {

PREDEPLOYMENT_CALIB_COEFFICIENT=ncvar_get(filenc_meta,"PREDEPLOYMENT_CALIB_COEFFICIENT")

PCC_CHLA=str_split(PREDEPLOYMENT_CALIB_COEFFICIENT[str_detect(PREDEPLOYMENT_CALIB_COEFFICIENT,"CHLA")],",")

SCALE=as.numeric(str_replace(PCC_CHLA[[1]][1],"SCALE_CHLA=",""))

DARK=as.numeric(str_replace(PCC_CHLA[[1]][2],"DARK_CHLA=",""))

return(list("DARK"=DARK,"SCALE"=SCALE))

} 
