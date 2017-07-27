#This function is used in the script 1pool_models.R
OnePool<-function(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch){
 
  Xt<-matrix(NA,length(data[,1]),2)
  
  r<-par[1]*(par[2]^(data$epiTemp-20))
  photo<-photo_switch*par[3]*exp(par[7]*data$PAR/70000)
  inpart<-exp(-par[6]*dens)
  orgpart<-exp(-par[6]*dens)
  if(dens_switch==0){
    inpart<-exp(0*dens)
    orgpart<-exp(0*dens)
  }
  
  if(nu_hDIC_switch==0){
    par[4]<-1
  }
  if(nu_hDOC_switch==0){
    par[5]<-1
  }
  
   if(run_type==1){ #this is the first time point index observations for DOC and CO2 were available in 2014
    s<-6
  }
  if(run_type==2){ #this is the first time point index observations for DOC and CO2 were available in 2015
    s<-10
  }
  
  for (t in s:length(data[,1])){
    if (t==s){
      Xt[t-1,1]<-data$dic[t-1]
      Xt[t-1,2]<-data$doc[t-1]
      
    }
    
    change_dic<-0
    change_doc<-0
    
    for (t2 in 1:1){
      change1<-(-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,1]-(data$kCO2[t-1]/data$thermo.depth[t-1])*Xt[t-1,1]+
                  r[t-1]*Xt[t-1,2]+((data$DICeq[t-1]*data$epiVol[t-1]/data$thermo.depth[t-1])*data$kCO2[t-1])+inpart[t-1]*data$dicIn[t-1]-
                  data$GPP[t-1]*(1-GPPrespired)+entr_switch*tohypo[t-1]*Xt[t-1,1]/data$epiVol[t-1]+entr_switch*par[4]*fromhypo[t-1]*data$hypo_dicInt[t-1]+photo[t-1]*Xt[t-1,2])/1
      
      change2<-(-r[t-1]*Xt[t-1,2]-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,2]+exudeTotal*data$GPP[t-1]+
                  orgpart[t-1]*data$docIn[t-1]+entr_switch*tohypo[t-1]*Xt[t-1,2]/data$epiVol[t-1]+entr_switch*par[5]*fromhypo[t-1]*data$hypo_docInt[t-1]-photo[t-1]*Xt[t-1,2])/1
      change_dic<-change_dic+change1
      change_doc<-change_doc+change2
    }
    
    Xt[t,1]<-Xt[t-1,1]+change_dic
    if(Xt[t,1]<0){
      Xt[t,1]<-0
    }
    Xt[t,2]<-Xt[t-1,2]+change_doc
  }
  return(Xt)
}