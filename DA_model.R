DA_model<-function(dat,p,D,t_est,n,col_nr,Rest,out,model){
  # Isolate element of interest
  M2<-dat[,c(1,col_nr)]
  M2<-M2[complete.cases(M2),]
  
  R_est<-Rest[1,col_nr] # Volumetric equilibrium coefficient for Zn into bioapatite (= Kp / p)
  
  # Find high values
  Mh<-M2[which(M2[,2] %in% tail(sort(M2[,2]),n)),]
  C1_est<-mean(Mh[,2])
  C1x_est<-mean(Mh[,1])
  # Find low values
  Ml<-M2[which(M2[,2] %in% head(sort(M2[,2]),n)),]
  C0_est<-mean(Ml[,2])
  C0x_est<-mean(Ml[,1])
  # Swap C0 and C1 if leaching is outward
  if(C1x_est>C0x_est){
    C0_est<-mean(Mh[,2])
    C1_est<-mean(Ml[,2])
  }
  # Find intermediate value
  Ci<-mean(c(C0_est,C1_est))
  Mi<-M2[which.min(abs(M2[,2]-Ci)),]
  # Estimate inflection point of diffusion front
  x0_est<-as.numeric(Mi[1])
  
  # Estimate pore fluid concentration C2 from outside concentration C1 corrected for C0 (negative in case of leachOUT)
  C2_est<-(C1_est-C0_est)/(R_est*p)
  # Define DA model that needs to be minimized. Concentration of C0 is subtracted from all data to add later to enable leachOUT simulations
  min.RSS<-function(par,data){
    with(data, sum(((par[2]-par[1])*erfc((data[,1]-par[3])/2*sqrt((R_est+1)/(D*par[4])))-(data[,2]-par[1]))^2))
  }
  # Create model based on first estimates, use opimalization because non-linear estimation does not work with error function
  par<-c(C0_est,C1_est-(C1_est-C0_est)/2,x0_est,t_est) # Bring C1 and C0 closer together for better initial values
  result<-optim(par,min.RSS,data=as.data.frame(M2),control=c(trace=F,maxit=10000),hessian=T)
  # Result with added C0 to simulate base concentrations
  model<-cbind(model,result$par[1]+(result$par[2]-result$par[1])*erfc((dat[,1]-result$par[3])/2*sqrt((R_est+1)/(D*result$par[4]))))
  colnames(model)[length(model[1,])]<-colnames(dat)[col_nr]
  # Plot result
  plot(dat[,c(1,col_nr)], main = colnames(dat)[col_nr],ylab="",xlab="");lines(model[,c(1,length(model[1,]))],col="red")
  # Calculate Pearson's correlation
  R2<-sqrt(1-result$value/sum((M2[,2]-mean(M2[,2]))^2))
  MM<-cbind(M2[,2],model[,2])
  PCorr<-cor.test(MM[,1],MM[,2],alternative="two.sided",method="pearson",exact=F)
  ppvalue<-PCorr$p.value
  pr<-as.numeric(PCorr$estimate)
  # Calculate Spearman correlation
  SpCorr<-cor.test(MM[,1],MM[,2],alternative="two.sided",method="spearman",exact=F)
  spvalue<-SpCorr$p.value
  rho<-as.numeric(SpCorr$estimate)
  # Calculate covariance matrix and standard errors of parameters
  Cov<-solve(result$hessian)
  SEpar<-sqrt(result$value/(length(M2[,1]-length(par)))*abs(diag(Cov)))
  # Calculate C2 (concentration of the fluid)
  C2<-result$par[2]/(p*R_est)
  C2_se<-SEpar[2]/(p*R_est)
  # Create output line
  if(C1x_est>C0x_est){
    INOUT<-"OUT"
  }else{INOUT<-"IN"}
  out<-rbind(out,c(colnames(dat)[col_nr],result$par,C2,SEpar,C2_se,pr,R2,ppvalue,rho,spvalue,INOUT))
  L<-list(out,model)
  return(L)
}

AppDAMod<-function(dat,p,D,t_est,n,Rest,sheet_name){
  # Create output file with headers
  model<-dat[,1]
  out<-c("parameter","C0","C1","x0","t","C2","C0_se","C1_se","x0_se","t_se","C2_se","pearsons r","pearsons R2","Pearsons p-value","Spearmans rho","Spearmans p-value","Direction of leaching")
  
  # Model all elements in dat and show model result fits
  dev.new()
  par(mfrow=c(5,5))
  for(i in 2:length(dat[1,])){
    try(nL<-DA_model(dat,p,D,t_est,n,i,Rest,out,model),TRUE)
    if(class(nL)=="try-error"){
      print(i)
      next
    }else{
      L<-nL
    }
    out<-as.data.frame(L[1],stringsAsFactors=F)
    model<-as.data.frame(L[2],stringsAsFactors=F)
  }
  write.xlsx(out,paste(sheet_name,"_model_output.xlsx"),col.names=F,row.names=F)
  write.xlsx(model,paste(sheet_name,"_model_data.xlsx"),col.names=T,row.names=F)
  L2<-list(out,model)
  return(L2)
}

# Define complementary error function
erfc <- function(x){2 * pnorm(x * sqrt(2), lower = FALSE)}

# Define DA model that needs to be minimized
min.RSS<-function(par,data){
  with(data, sum((par[1]+p*R*par[2]*erfc((data[,1]-par[3])/2*sqrt((R+1)/(D*par[4])))-data[,2])^2))
}

# Define point-moving-average function
pma<-function(x,i,n){
  ts<-cbind(x[,1],x[,i])
  ts<-cbind(ts,rep(NA,length(x[,1])))
  for (a in ((n/2+0.5):(length(x[,1])+0.5-(n/2)))) {
    ts[a,3] = mean(ts[(a-(n/2)+1):(a+(n/2)),2],na.rm=T)
  }		
  return(ts)
}

DA_model_MC_loop<-function(N,data,p,D,t_est,n,Rest,name){
  dev.new()
  par(mfrow=c(5,5))
  for(i in seq(2,length(data)-1,2)){
    stats<-DA_model_MC(N,data,p,D,t_est,n,i,Rest,F)
    model_stat<-as.data.frame(stats[[1]])
    write.xlsx(model_stat,paste(name,"_model_stat",colnames(data)[i],".xlsx"))
    par_stat<-as.data.frame(stats[[2]])
    write.xlsx(par_stat,paste(name,"_par_stat",colnames(data)[i],".xlsx"))
  }
}

# Monte Carlo simulation of model skipping errors
DA_model_MC<-function(N,data,p,D,t_est,n,col_nr,Rest,export_full){
  plot(data[,c(1,col_nr)], main = colnames(data)[col_nr],ylab="",xlab="")
  x<-data[,1]
  y<-data[,col_nr]
  y_sd<-data[,col_nr+1]
  model<-x
  parls<-c("C0","C1","x0","t")
  MP<-list(model,parls)
  # Subsample new Y-values
  y_mod<-vector()
  for(i in 1:length(y)){
    y_mod<-as.data.frame(rbind(y_mod,rnorm(N,y[i],y_sd[i])))
  }
  y_mod<-cbind(x,y_mod)
  colnames(y_mod)<-c("x",seq(1,N,1))
  R_est<-Rest[1,col_nr/2+1]
  for(i in 1:N){
    tryCatch({
      MP<-DA_model2(y_mod,p,D,t_est,n,i+1,R_est,model,parls)
    },error=function(e){
      cat("Error:",conditionMessage(e),"\n");
      MP<-list(model,parls)
    })
    model<-as.data.frame(MP[[1]])
    parls<-as.data.frame(MP[[2]])
  }
  # Calculate stats
  model_mean<-apply(model,1,mean)
  model_sd<-apply(model,1,sd)
  model_c95_up<-model_mean+qnorm(0.975)*model_sd
  model_c95_down<-model_mean-qnorm(0.975)*model_sd
  model_c95_up_se<-model_mean+qnorm(0.975)*model_sd/sqrt(N)
  model_c95_down_se<-model_mean-qnorm(0.975)*model_sd/sqrt(N)
  model_stat<-cbind(x,model_mean,model_sd,model_c95_up,model_c95_down, model_c95_up_se, model_c95_down_se)
  colnames(model_stat)<-c("X","Modelled mean value","Modelled standard deviation","95% cl up","95% cl down","95% cl up SE","95% cl down SE")
  
  # Remove factors from parls
  rownames(parls)<-parls[,1]
  parls<-parls[,-1]
  parls<-data.frame(lapply(parls, as.character), stringsAsFactors=FALSE)
  parls<-data.frame(lapply(parls, as.numeric))
  
  par_mean<-apply(parls[,-1],1,mean)
  par_sd<-apply(parls[,-1],1,sd)
  par_se<-par_sd/sqrt(N)
  par_stat<-cbind(c("C0","C1","x0","t"),par_mean,par_sd,par_se)
  
  #plot(data[,c(1,col_nr)], main = colnames(data)[col_nr],ylab="",xlab="")
  lines(model_stat[,c(1,2)])
  lines(model_stat[,c(1,4)],col="red",lty=3)
  lines(model_stat[,c(1,5)],col="red",lty=3)
  lines(model_stat[,c(1,6)],col="green",lty=3)
  lines(model_stat[,c(1,7)],col="green",lty=3)
  stats<-list(model_stat,par_stat)
  if(export_full){
    full_results<-list(model_stat,par_stat,y_mod,model,parls)
    return(full_results)
  }else{
    return(stats)
  }
}

# Model best fitting DA curve
DA_model2<-function(dat,p,D,t_est,n,col_nr,R_est,model,parls){
  # Isolate element of interest
  M2<-dat[,c(1,col_nr)]
  M2<-M2[complete.cases(M2),]
  
  # Find high values
  Mh<-M2[which(M2[,2] %in% tail(sort(M2[,2]),n)),]
  C1_est<-mean(Mh[,2])
  C1x_est<-mean(Mh[,1])
  # Find low values
  Ml<-M2[which(M2[,2] %in% head(sort(M2[,2]),n)),]
  C0_est<-mean(Ml[,2])
  C0x_est<-mean(Ml[,1])
  # Swap C0 and C1 if leaching is outward
  if(C1x_est>C0x_est){
    C0_est<-mean(Mh[,2])
    C1_est<-mean(Ml[,2])
  }
  # Find intermediate value
  Ci<-mean(c(C0_est,C1_est))
  Mi<-M2[which.min(abs(M2[,2]-Ci)),]
  # Estimate inflection point of diffusion front
  x0_est<-as.numeric(Mi[1])
  
  # Estimate pore fluid concentration C2 from outside concentration C1 corrected for C0 (negative in case of leachOUT)
  C2_est<-(C1_est-C0_est)/(R_est*p)
  # Define DA model that needs to be minimized. Concentration of C0 is subtracted from all data to add later to enable leachOUT simulations
  min.RSS<-function(par,data){
    with(data, sum(((par[2]-par[1])*erfc((data[,1]-par[3])/2*sqrt((R_est+1)/(D*par[4])))-(data[,2]-par[1]))^2))
  }
  # Create model based on first estimates, use opimalization because non-linear estimation does not work with error function
  par<-c(C0_est,C1_est-(C1_est-C0_est)/2,x0_est,t_est) # Bring C1 and C0 closer together for better initial values
  result<-optim(par,min.RSS,data=as.data.frame(M2),control=c(trace=F,maxit=10000),hessian=T)
  # Result with added C0 to simulate base concentrations
  model<-cbind(model,result$par[1]+(result$par[2]-result$par[1])*erfc((dat[,1]-result$par[3])/2*sqrt((R_est+1)/(D*result$par[4]))))
  colnames(model)[length(model[1,])]<-colnames(dat)[col_nr]
  # Save modelled parameters
  parls<-cbind(parls,result$par)
  colnames(parls)[length(parls[1,])]<-colnames(dat)[col_nr]
  # Plot result
  #lines(model[,c(1,length(model[1,]))],col="red")
  MP<-list(model,parls)
  return(MP)
}
