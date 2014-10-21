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

######
# estimate depth distribution of seagrass
max_est <- function(dat_in, depth_var = 'depth', sg_var = 'seagrass', 
  dat_out = F){
  
	# order by depth, assumes column is negative
  dat_in <- dat_in[order(dat_in[, depth_var], decreasing = T), ]
	dat_in$depth <- dat_in[, depth_var]
	
	# cumulative sum of pts with all seagrass and all points
	# assumes NA is empty
	sg_pts <- rev(table(dat_in[!is.na(dat_in[, sg_var]), depth_var]))
	sg_pts <- data.frame(Depth = names(sg_pts), sg_pts = sg_pts,
		sg_cum = cumsum(sg_pts), row.names = 1:length(sg_pts))
	dep_pts <- rev(table(dat_in[, depth_var]))
	dep_pts <- data.frame(Depth = names(dep_pts), dep_pts = dep_pts, 
		dep_cum = cumsum(dep_pts), row.names = 1:length(dep_pts))
	
	# combine all pts and seagrass pts, depth as numeric
	pts <- merge(dep_pts, sg_pts, by = 'Depth', all.x = T)
	pts$Depth <- as.numeric(as.character(pts$Depth))
	pts$sg_prp <- with(pts, sg_pts/dep_pts)
	
	# add slope ests to pts
	pts$dep_slo <- with(pts, c(NA, -1 * diff(dep_cum)/diff(Depth)))
	pts$sg_slo <- with(pts, c(NA, -1 * diff(sg_cum)/diff(Depth)))
	
	# return cumulative data if T
	if(dat_out) return(pts)
	
	# get max (95th) and median (50th) depth for all seagrass
	max_depth_all <- dat_in[which.min(abs(dat_in[, 'all'] - 0.95)), 'depth']
	med_depth_all <- dat_in[which.min(abs(dat_in[, 'all'] - 0.5)), 'depth']
	
	# get max (95th) and median (50th) depth for continuous seagrass
	max_depth_cont <- dat_in[which.min(abs(dat_in[, 'cont'] - 0.95)), 'depth']
	med_depth_cont <- dat_in[which.min(abs(dat_in[, 'cont'] - 0.5)), 'depth']
		
  # return output
  out <- c(zmax_all = max_depth_all, z50_all = med_depth_all, 
  	zmax_cont = max_depth_cont, z50_cont = med_depth_cont)
  return(out)
  
}

######
# functions for seagrass depth of col estimates

######
# function for estimating depth of colonization
# also used for plots
# 'dat_in' is data from 'buff_ext'
# 'depth_var' is name of depth column in input data
# 'sg_var' is name of seagrass column in input data
# 'thresh' is numeric threshold value for estimating depth of col
doc_est <- function(dat_in, depth_var = 'Depth', sg_var = 'Seagrass',
  thresh = 0.1){
  
  # order by depth, assumes column is negative
  dat_in <- dat_in[order(dat_in[, depth_var], decreasing = T), ]
	dat_in$depth <- dat_in[, depth_var]
	
	# cumulative sum of pts with all seagrass and all points
	# assumes NA is empty
	sg_pts <- table(dat_in[!is.na(dat_in[, sg_var]), depth_var])
	sg_pts <- data.frame(Depth = names(sg_pts), sg_pts = as.numeric(sg_pts),
		sg_cum = cumsum(sg_pts), row.names = 1:length(sg_pts))
	dep_pts <- table(dat_in[, depth_var])
	dep_pts <- data.frame(Depth = names(dep_pts), dep_pts = as.numeric(dep_pts), 
		dep_cum = cumsum(dep_pts), row.names = 1:length(dep_pts))
	
	# combine all pts and seagrass pts, depth as numeric
	pts <- merge(dep_pts, sg_pts, by = 'Depth', all.x = T)
	pts$Depth <- as.numeric(as.character(pts$Depth))
	# output
  pts$sg_prp <- with(pts, sg_pts/dep_pts)
	
	##
	# estimate a logistic growth function for the data
	resps <- c('dep_cum', 'sg_cum')
	pred_ls <- vector('list', length(resps))
	names(pred_ls) <- resps
	for(resp in resps){
		
		##
		# estimates of starting parameters
		
		# logistic growth
		Asym <- max(pts[, resp], na.rm = T)
		xmid <- median(pts$Depth, na.rm = T)
		scal <- quantile(pts$Depth, 0.75, na.rm = T) - xmid
		form_in <- substitute(x ~ SSlogis(Depth, Asym,  xmid, scal), list(x = as.name(resp)))

#			# Gompertz
# 		Asym <- max(pts[, resp], na.rm = T)
# 		b2 <- median(pts$Depth, na.rm = T)
# 		b3 <- b3
# 		form_in <- substitute(x ~ SSgompertz(Depth, Asym, b2, b3), list(x = as.name(resp)))

		# model
		mod <- try(nls(form_in, data = pts, na.action = na.exclude))
		
    # values for prediction
		dep_rng <- range(pts[, 'Depth'], na.rm = T)
		new.x <- seq(dep_rng[1], dep_rng[2], length = 100)
	
    # return NAs if model fail, else get model predictions
    if('try-error' %in% class(mod)) {
      pred_ls[[resp]] <- rep(NA_real_, length = length(new.x))
    } else {
		pred_ls[[resp]] <- as.numeric(predict(mod, 
			newdata = data.frame(Depth = new.x)))
    }
		
	}
	
	preds <- data.frame(Depth = new.x, do.call('cbind', pred_ls))

	# add slope ests to pts, use differences
	preds$dep_slo <- with(preds, c(NA, diff(dep_cum)/diff(Depth)))
	preds$sg_slo <- with(preds, c(NA, diff(sg_cum)/diff(Depth)))
	
	# add threshold data based on proportion of dep_slo 
	threshs <- sapply(1:length(thresh), 
		FUN = function(x) thresh[x] * preds[, 'dep_slo']
		)
	threshs <- data.frame(threshs)
	names(threshs) <- paste('Threshold', thresh)
	
  # output
	preds <- data.frame(preds, threshs)
	
	# calculate depth of col
	doc <- sapply(thresh, 
		FUN = function(x){
			
			col <- preds[, grep(x, names(preds))]
			ind <- with(preds,  sg_slo <= col)
      ind <- max(which(!ind)) + 1
			preds[ind, 'Depth']
			
			}
		)
  # output
	names(doc) <- thresh

  # all output
	return(list(data = pts, preds = preds, ests = doc))
	  
}

#######
# function for creating random grid of points, bounded by polygon extent
# taken from ibi sampling manuscript functions
# 'clip_poly' is shapefile input object
# 'spacing' is spacing between points, as degrees
grid_est <- function(clip_poly, spacing = 0.03){
  
  if(!'SpatialPolygonsDataFrame' %in% class(clip_poly))
    stop('clip_poly must be of class SpatialPolygonsDataFrame')
  
  library(sp) 
  
  # extent of shapefile
  extent <- summary(clip_poly)$bbox
  
  # buffer of shapefile and random starting x/y locs
  add.on <- apply(extent, 1, diff) * 0.3
  rand <- runif(2, 0, spacing)
  
  # random points within rectangular extent
  pts<-{
    x.vals<-seq(extent[1, 1] - add.on['x'], extent[1, 2] + add.on['x'], by = spacing) + rand[1]
    y.vals<-seq(extent[2, 1] - add.on['y'], extent[2, 2] + add.on['y'], by = spacing) + rand[2]
    expand.grid(x.vals, y.vals)
  }
  
  # clip by clip_poly and return
  sel <- !is.na(SpatialPoints(pts) %over% clip_poly)[, 1]
  
  return(SpatialPoints(pts)[sel, ])
  
}

######
# function extracts bathymetric seagrass pts withing a distance from a pt
# 'pts' is spatial points to extract
# 'center' is pt from which buffer extends
# 'buff' is radius of buffer in dec degrees
buff_ext <- function(pts, center, buff = 0.03){
  
  # sanity checks
  if(!any(c('SpatialPointsDataFrame', 'SpatialPoints') %in% class(pts)))
    stop('pts must be of class SpatialPointsDataFrame or SpatialPoints')
  
  if(!any(c('SpatialPointsDataFrame', 'SpatialPoints') %in% class(center)))
    stop('center must be of class SpatialPointsDataFrame or SpatialPoints')
  
  library(rgeos)
  library(sp)
  
  # create buffer
  buffer <- gBuffer(center, width = buff)
  
  # index of pts in buffer
  sel <- !is.na(pts %over% buffer)
  
  if(sum(sel) == 0) stop('No points in buffer')
  
  return(pts[sel, ])
  
}

######
# function for plotting repeating polygons in ggplot
# 'flag.in' is vector indicating binomial variable for polys
# 'flag.in' can be factor or numeric (w/ two values)
# 'dat' is data.frame with 'flag.in'
# 'for_leg' is used to manually change color for polygons, will add legend
# output is geom object
poly.fun<-function(flag.in,dat, for_leg = F){

  require(reshape2)
  require(ggplot2)
  
  #for flag bias
  if(class(flag.in) == 'numeric'){ 
    neg.dates<-with(
      dat,
      DateTimeStamp[which(flag.in < 0)]
      )
    tz<-attr(neg.dates,'tzone')
    diffs<-c(0,diff(as.numeric(neg.dates)))
    strt.ind<-c(1,which(diffs>1800))
    end.ind<-c(strt.ind-1,length(flag.in))
    comb<-paste(neg.dates[strt.ind],neg.dates[end.ind],sep="\t")
    
    if(grepl('NA',comb[length(comb)])) 
      comb[length(comb)]<-gsub('NA',neg.dates[length(neg.dates)],comb[length(comb)])
  
    comb<-do.call('rbind',strsplit(comb,'\t'))
    comb<-cbind(comb,comb[,2],comb[,1])
  
    x.vals<-suppressMessages(melt(sapply(1:nrow(comb), 
      function(x) comb[x,],simplify=F))$value)
    x.vals<-as.POSIXct(as.character(x.vals),tz,
      format='%Y-%m-%d %H:%M:%S')
    y.vals<-rep(c(-1000,-1000,1000,1000),length=length(x.vals)) 
    Antagonistic<-rep(1:(length(x.vals)/4),each=4)
    polys<-data.frame(x.vals,y.vals,grp=Antagonistic)
    
    }

  #for sunset/rise
  if(class(flag.in) == 'factor'){
    
    plo.dates<-unique(dat[,c('solar','value')])
    
    if(plo.dates$solar[1] == 'sunset') 
      plo.dates<-plo.dates[-1,]
    if(plo.dates$solar[nrow(plo.dates)] == 'sunrise')
      plo.dates<-rbind(
        plo.dates,
        data.frame(solar='sunset',value=max(dat$DateTimeStamp))
        )
    
    plo.dates$inds<-rep(1:(nrow(plo.dates)/2),each=2)
    tz<-attr(plo.dates$value,'tzone')
    plo.dates$value <- as.character(plo.dates$value)
    plo.dates<-dcast(plo.dates,inds~solar,value.var='value')
    plo.dates<-with(plo.dates,
      data.frame(sunrise,sunset,sunset,sunrise)
      )
   
    x.vals <- sapply(1:nrow(plo.dates), 
           function(x) plo.dates[x,],simplify=F)
    x.vals<-suppressMessages(
      melt(x.vals, measure.vars = names(x.vals[[1]]))$value
      )
    x.vals<-as.POSIXct(x.vals, tz, origin = '1970-01-01')
    y.vals<-rep(c(-1000,-1000,1000,1000),nrow(plo.dates))
    Day<-as.character(trunc(x.vals,'days'))
    polys<-data.frame(x.vals,y.vals,grp=Day)

    }
  
  if(for_leg){
    out<-geom_polygon(data=polys,aes(x.vals,y.vals,group=grp, fill = 'grp'), 
      alpha=0.6)
  } else {
   out<-geom_polygon(data=polys,aes(x.vals,y.vals,group=grp), fill = "#EBCC2A",
      alpha=0.6)
  }
    
  
  return(out)
  
  }

#functions for NEM processing of NERRS data
#created Dec. 2013 by M. Beck, adapted from 'spam_NEM_fun.r' and M. Murrell

#funcion that splits dataset into 24hr days based on sunrise
#merge with original data
met.day.fun<-function(dat.in, stat.in, 
  meta.path = 'M:/wq_models/SWMP/sampling_stations.csv'
  ){

  require(StreamMetabolism)  #for sunrise.set function
  
  if(!exists('dat.meta')) 
    dat.meta<-read.csv(meta.path,header=T)
  
  stat.meta<-toupper(paste0(stat.in,'WQ'))
  stat.meta<-dat.meta[grep(stat.meta,toupper(dat.meta$Station.Code)),]
  
  # all times are standard - no DST!
  gmt.tab<-data.frame(
    gmt.off=c(-4,-5,-6,-8,-9),
    tz=c('America/Virgin', 'America/Jamaica', 'America/Regina',
      'Pacific/Pitcairn', 'Pacific/Gambier'),
    stringsAsFactors=F
    )
  
  #get sunrise/sunset times using sunrise.set function from StreamMetabolism
  lat<-stat.meta$Latitude
  long<--1*stat.meta$Longitude
  GMT.Offset<-stat.meta$GMT.Offset
  tz<-gmt.tab[gmt.tab$gmt.off==GMT.Offset,'tz']
  start.day<-format(dat.in$DateTimeStamp[which.min(dat.in$DateTimeStamp)]-(60*60*24),format='%Y/%m/%d')
  tot.days<-1+length(unique(as.Date(dat.in$DateTimeStamp)))
  
  #ss.dat is matrix of sunrise/set times for each days  within period of obs
  ss.dat<-suppressWarnings(sunrise.set(lat,long,start.day,tz,tot.days))
  
  #remove duplicates, sometimes sunrise.set screws up
  ss.dat<-ss.dat[!duplicated(strftime(ss.dat[,1],format='%Y-%m_%d')),]
  ss.dat<-data.frame(
    ss.dat,
    met.date=as.Date(ss.dat$sunrise,tz=tz)
    )
  ss.dat<-melt(ss.dat,id.vars='met.date')
  if(!"POSIXct" %in% class(ss.dat$value))
    ss.dat$value<-as.POSIXct(ss.dat$value, origin='1970-01-01',tz=tz)
  ss.dat<-ss.dat[order(ss.dat$value),]
  ss.dat$day.hrs<-unlist(lapply(
    split(ss.dat,ss.dat$met.date),
    function(x) rep(as.numeric(x[2,'value']-x[1,'value']),2) 
    ))
  
  #matches is vector of row numbers indicating starting value that each
  #unique DateTimeStamp is within in ss.dat
  #output is meteorological day matches appended to dat.in
  matches<-findInterval(dat.in$DateTimeStamp,ss.dat$value)
  data.frame(dat.in,ss.dat[matches,])
      
  }

######
# get legend from an existing ggplot object
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
