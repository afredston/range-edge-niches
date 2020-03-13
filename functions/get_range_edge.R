# this function is adapted from plot_range_edge but will not make plots for speed and ease of storage of many species results 
# it requires reshape2

# Sdreport=fit$parameter_estimates$SD
# Report=fit$Report
# TmbData=fit$data_list
# Obj=fit$tmb_list$Obj
# working_dir=paste0(getwd(),"/")
# Year_Set=fit$year_labels
# n_samples=100
# quantiles=c(0.1,0.5,0.9) 

get_range_edge = function( fit.model, 
                           strata_names=NULL,
                            category_names=NULL,
                           working_dir=paste0(getwd(),"/"), 
                           quantiles=c(0.05,0.95), 
                           n_samples=100,
                            interval_width=1, 
                           width=NULL, 
                           height=NULL,
                           calculate_relative_to_average=FALSE,
                           ...){
  
  # Unpack
  
  Sdreport=fit.model$parameter_estimates$SD 
  Report=fit.model$Report 
  TmbData=fit.model$data_list 
  Obj=fit.model$tmb_list$Obj
  year_labels=fit$year_labels
  
  # # Informative errors
  # if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")
  # if( !("jointPrecision" %in% names(Sdreport))) stop("jointPrecision not present in Sdreport; please re-run with `getJointPrecision=TRUE`")
  # if( any(quantiles<0) | any(quantiles>1) ) stop("Please provide `quantiles` between zero and one")
  # if( all(TmbData$Z_gm==0) ) stop("Please re-run with 'Options['Calculate_Range']=TRUE' to calculate range edges")
  # 
  # # Which parameters
  # if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
  #   # SpatialDeltaGLMM
  #   stop("Not implemente")
  # }
  # if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
  #   # VAST Version < 2.0.0
  #   stop("Not implemente")
  # }
  # if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
  #   # VAST Version >= 2.0.0
  #   TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  # }
  
  # Default inputs
#  if( is.null(Year_Set)) Year_Set = 1:TmbData$n_t
#  if( is.null(Years2Include) ) Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c
  if( is.null(colnames(TmbData$Z_gm)) ){
    m_labels = paste0("axis",1:ncol(TmbData$Z_gm))
  }else{
    m_labels = colnames(TmbData$Z_gm)
  }
  
  
  ##### Local function
  D_gcyr = sample_variable( Sdreport=Sdreport, Obj=Obj, variable_name="D_gcy", n_samples=n_samples )
  
  # Calculate quantiles from observed and sampled densities D_gcy
  E_zctm = array(NA, dim=c(length(quantiles),dim(Report$D_gcy)[2:3],ncol(TmbData$Z_gm)) )
  E_zctmr = array(NA, dim=c(length(quantiles),dim(Report$D_gcy)[2:3],ncol(TmbData$Z_gm),n_samples) )
  Mean_cmr = array(NA, dim=c(dim(Report$D_gcy)[2],ncol(TmbData$Z_gm),n_samples) )
  prop_zctm = array(NA, dim=c(dim(Report$D_gcy)[1:3],ncol(TmbData$Z_gm)) )
  prop_zctmr = array(NA, dim=c(dim(Report$D_gcy)[1:3],ncol(TmbData$Z_gm),n_samples) )
  for( rI in 0:n_samples ){
    for( mI in 1:ncol(TmbData$Z_gm) ){
      order_g = order(TmbData$Z_gm[,mI], decreasing=FALSE)
      if(rI==0) prop_zctm[,,,mI] = apply( Report$D_gcy, MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
      if(rI>=0) prop_zctmr[,,,mI,rI] = apply( D_gcyr[,,,rI,drop=FALSE], MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
      
      # Calculate edge
      for( cI in 1:dim(E_zctm)[2] ){
        if(rI>=1){
          if( calculate_relative_to_average==TRUE ){
            Mean_cmr[cI,mI,rI] = weighted.mean( as.vector(TmbData$Z_gm[,mI]%o%rep(1,dim(Report$D_gcy)[3])), w=as.vector(D_gcyr[,cI,,rI]) )
          }else{
            Mean_cmr[cI,mI,rI] = 0
          }
        }
        for( zI in 1:dim(E_zctm)[1] ){
          for( tI in 1:dim(E_zctm)[3] ){
            if(rI==0){
              index_tmp = which.min( (prop_zctm[,cI,tI,mI]-quantiles[zI])^2 )
              E_zctm[zI,cI,tI,mI] = TmbData$Z_gm[order_g[index_tmp],mI]
            }
            if(rI>=1){
              index_tmp = which.min( (prop_zctmr[,cI,tI,mI,rI]-quantiles[zI])^2 )
              E_zctmr[zI,cI,tI,mI,rI] = TmbData$Z_gm[order_g[index_tmp],mI] - Mean_cmr[cI,mI,rI]
            }
          }}
      }
    }}
  SE_zctm = apply( E_zctmr, MARGIN=1:4, FUN=sd )
  Edge_zctm = abind::abind( "Estimate"=E_zctm, "Std. Error"=SE_zctm, along=5 )
  dimnames(Edge_zctm)[[1]] = paste0("quantile_",quantiles)
  
  # transform matrix into a dataframe 
  
  Edge_df <- reshape2::melt(Edge_zctm)
  Edge_df$Var2 <- NULL
  
  # replace axis numbers with names 
  for(i in 1:length(Z_gm_axes)) {
    Edge_df$Var4 <- gsub(paste0(i), paste0(Z_gm_axes[i]), Edge_df$Var4)
  }

  colnames(Edge_df) <- c("quantile","year","axis","quantity","value")
  Edge_df$year <- year_labels[Edge_df$year]  # DANGER--would prefer to carry through real year values
  
  # eliminate unwanted years from Edge_df
  
  Edge_df <- Edge_df[Edge_df$year %in% Years2Include,]
  
  # Return list of stuff
  
  return(Edge_df)
}
