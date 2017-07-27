#This function is used in the script 2pool_models.R
TwoPool<-function(data,par,run_type,entr_switch,nu_hDIC_switch,nu_hDOC_switch,dens_switch,photo_switch,rDOC_to_lDOC_switch,var_iDOC_lability){
 
  Xt<-matrix(NA,length(data[,1]),3)
  
  if(nu_hDIC_switch==0){
    par[9]<-1 
  }
  if(nu_hDOC_switch==0){
    par[10]<-1
  }
  
   if(run_type==1){ #this is the first time point index observations for DOC and CO2 were available in 2014
    s<-6
  }
  if(run_type==2){ #this is the first time point index observations for DOC and CO2 were available in 2015
    s<-10
  }
  
  r<-par[1]*(par[2]^(data$epiTemp-20))
  r2<-par[3]*(par[2]^(data$epiTemp-20))
  inpart<-exp(-par[8]*dens)
  orgpart<-exp(-par[8]*dens)
  
  if(dens_switch==0){
    inpart<-exp(0*dens)
    orgpart<-exp(0*dens)
  }
 
  if (rDOC_to_lDOC_switch==0){
    par[7]<-0
  } else {
    par[7]<-par[7]
  }


  photoL<-photo_switch*par[13]*exp(par[15]*data$PAR/70000)
  photoR<-photo_switch*par[14]*exp(par[16]*data$PAR/70000)
  
  if (var_iDOC_lability==0){
  for (t in s:length(data[,1])){
    if (t==s){
      Xt[t-1,1]<-data$dic[t-1]
      Xt[t-1,2]<-par[4]*data$doc[t-1]
      Xt[t-1,3]<-(1-par[4])*data$doc[t-1]
    }
   
    change_dic<-0
    change_doc<-0
    change_docr<-0
    
    for (t2 in 1:1){
      change1<-(-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,1]-(data$kCO2[t-1]/data$thermo.depth[t-1])*Xt[t-1,1]+
                  r[t-1]*Xt[t-1,2]+(1-par[7])*r2[t-1]*Xt[t-1,3]+((data$DICeq[t-1]*data$epiVol[t-1]/data$thermo.depth[t-1])*data$kCO2[t-1])+inpart[t-1]*data$dicIn[t-1]-
                  data$GPP[t-1]*(1-GPPrespired)+entr_switch*tohypo[t-1]*Xt[t-1,1]/data$epiVol[t-1]+entr_switch*par[9]*fromhypo[t-1]*data$hypo_dicInt[t-1]+photoL[t-1]*Xt[t-1,2])/1
      
      
      change2<-(-r[t-1]*Xt[t-1,2]-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,2]+exudeTotal*data$GPP[t-1]+
                  (1-par[5])*orgpart[t-1]*data$docIn[t-1]+entr_switch*tohypo[t-1]*Xt[t-1,2]/data$epiVol[t-1]+entr_switch*par[10]*par[6]*fromhypo[t-1]*data$hypo_docInt[t-1] -
                  photoL[t-1]*Xt[t-1,2]+photoR[t-1]*Xt[t-1,3]+par[7]*r2[t-1]*Xt[t-1,3])/1
      
      change3<-(-r2[t-1]*Xt[t-1,3]-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,3]+entr_switch*tohypo[t-1]*Xt[t-1,3]/data$epiVol[t-1]+entr_switch*par[10]*(1-par[6])*fromhypo[t-1]*data$hypo_docInt[t-1]+
                  par[5]*orgpart[t-1]*data$docIn[t-1]-photoR[t-1]*Xt[t-1,3])/1
      
      
      change_dic<-change_dic+change1
      change_doc<-change_doc+change2
      change_docr<-change_docr+change3
    }
    
    Xt[t,1]<-Xt[t-1,1]+change_dic
    Xt[t,2]<-Xt[t-1,2]+change_doc
    Xt[t,3]<-Xt[t-1,3]+change_docr
  }
  }
  
  if (var_iDOC_lability==1){
    for (t in s:length(data[,1])){
      if (t==s){
        Xt[t-1,1]<-data$dic[t-1]
        Xt[t-1,2]<-par[4]*data$doc[t-1]
        Xt[t-1,3]<-(1-par[4])*data$doc[t-1]
      }
      
      change_dic<-0
      change_doc<-0
      change_docr<-0
      if(t<=par[11]){
        par[5]<-par[5]
      }
      if(t>day){
        par[5]<-par[12]
      }
      for (t2 in 1:1){
        change1<-(-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,1]-(data$kCO2[t-1]/data$thermo.depth[t-1])*Xt[t-1,1]+
                    r[t-1]*Xt[t-1,2]+(1-par[7])*r2[t-1]*Xt[t-1,3]+((data$DICeq[t-1]*data$epiVol[t-1]/data$thermo.depth[t-1])*data$kCO2[t-1])+inpart[t-1]*data$dicIn[t-1]-
                    data$GPP[t-1]*(1-GPPrespired)+entr_switch*tohypo[t-1]*Xt[t-1,1]/data$epiVol[t-1]+entr_switch*par[9]*fromhypo[t-1]*data$hypo_dicInt[t-1]+photoL[t-1]*Xt[t-1,2])/1
        
        
        change2<-(-r[t-1]*Xt[t-1,2]-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,2]+exudeTotal*data$GPP[t-1]+
                    (1-par[5])*orgpart[t-1]*data$docIn[t-1]+entr_switch*tohypo[t-1]*Xt[t-1,2]/data$epiVol[t-1]+entr_switch*par[10]*par[6]*fromhypo[t-1]*data$hypo_docInt[t-1] -
                    photoL[t-1]*Xt[t-1,2]+photoR[t-1]*Xt[t-1,3]+par[7]*r2[t-1]*Xt[t-1,3])/1
        
        change3<-(-r2[t-1]*Xt[t-1,3]-(data$QoutInt[t-1]/data$epiVol[t-1])*Xt[t-1,3]+entr_switch*tohypo[t-1]*Xt[t-1,3]/data$epiVol[t-1]+entr_switch*par[10]*(1-par[6])*fromhypo[t-1]*data$hypo_docInt[t-1]+
                    par[5]*orgpart[t-1]*data$docIn[t-1]-photoR[t-1]*Xt[t-1,3])/1
        
        
        change_dic<-change_dic+change1
        change_doc<-change_doc+change2
        change_docr<-change_docr+change3
      }
      
      Xt[t,1]<-Xt[t-1,1]+change_dic
      Xt[t,2]<-Xt[t-1,2]+change_doc
      Xt[t,3]<-Xt[t-1,3]+change_docr
    }
  }
  Xt<-cbind(Xt[,1],Xt[,2]+Xt[,3])
  return(Xt)
}