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