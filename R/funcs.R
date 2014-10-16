######
#function for getting regression weights
#'wt.vars' is name of three variables to weight
#'ref.in' is row of dat.in that is used as reference
#'dat.in' is data to get weights from
#'wins' are the windows for the three wt.vars
#'all' will return all weights, rather than the product of all three
#'min.obs' uses window widening if less than 100 non-zero weights are found
wt.fun<-function(ref.in,dat.in,
  wt.vars=c('month.num','year','sal.ref'),
  wins=list(0.5,10,NULL),
  all=F,
  min.obs=T){
  
  #sanity check
  if(sum(wt.vars %in% names(dat.in)) != length(wt.vars))
    stop('Weighting variables must be named in "dat.in"')
  
  #windows for each of three variables
  wins_1<-wins[[1]]
  wins_2<-wins[[2]]
  wins_3<-wins[[3]]
  
  if(is.null(wins[[3]])) wins_3<-diff(range(dat.in[,wt.vars[3]]))/2
  
  #weighting tri-cube function
  wt.fun.sub<-function(dat.cal,ref,win,mo=F){
    
    dist.val<-abs(dat.cal-ref)
    
    if(mo){
      dist.val<-pmin(
        ref+1-dat.cal,
        dat.cal+1-ref,
        dist.val
        )
      }
    
    if(dist.val<=win) return((1-(dist.val/win)^3)^3)
    
    else return(0)
      
    }

  #reference (starting) data
  ref_1<-as.numeric(ref.in[,wt.vars[1]])
  ref_2<-as.numeric(ref.in[,wt.vars[2]])
  ref_3<-as.numeric(ref.in[,wt.vars[3]])

  #weights for each observation in relation to reference
  wts_1<-sapply(as.numeric(dat.in[,wt.vars[1]]),wt.fun.sub,ref=ref_1,win=wins_1,mo=T)
  wts_2<-sapply(as.numeric(dat.in[,wt.vars[2]]),wt.fun.sub,ref=ref_2,win=wins_2)
  wts_3<-sapply(as.numeric(dat.in[,wt.vars[3]]),wt.fun.sub,ref=ref_3,win=wins_3)
  out<-wts_1*wts_2*wts_3
  
  gr.zero<-sum(out>0)
  #cat('   Number of weights greater than zero =',gr.zero,'\n')
  
  while(gr.zero<100){
    
    wins_1<-0.1*wins_1 + wins_1
    wins_2<-0.1*wins_2 + wins_2
    wins_3<-0.1*wins_3 + wins_3
    
    #weights for each observation in relation to reference
    wts_1<-sapply(as.numeric(dat.in[,wt.vars[1]]),wt.fun.sub,ref=ref_1,win=wins_1,mo=T)
    wts_2<-sapply(as.numeric(dat.in[,wt.vars[2]]),wt.fun.sub,ref=ref_2,win=wins_2)
    wts_3<-sapply(as.numeric(dat.in[,wt.vars[3]]),wt.fun.sub,ref=ref_3,win=wins_3)
    
    out<-wts_1*wts_2*wts_3
    
    gr.zero<-sum(out>0)
    
    #cat('   Number of weights greater than zero',gr.zero,'\n')
    
    }
  
  #return all weights if T
  if(all){
    out<-data.frame(wts_1,wts_2,wts_3)
    names(out)<-wt.vars
    return(out)
    }
  
  #final weights are product of all three
  out
  
  }
