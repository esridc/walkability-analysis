bind_up_indicators = function(AOI){
  s = Sys.time()
  folders = list.dirs(paste("/Users/aaro9563/Documents/01 Walkability/", AOI, "/DUC Input", sep=""))[2:5]
  cats = c("grocery","health","schools","transit")
  shapes = list()
  for(i in 1:4){
    shapes[[cats[i]]] = list.files(folders[i], pattern=".shp$", full.names=T)
  }
  rm(i, folders)
  sfs = lapply(shapes, function(x){
    y = read_sf(x)
    y = st_transform(y, 102005)
    if(any(is.na(st_dimension(y)))){
      y = y[-which(is.na(st_dimension(y))),]
    }
    y = y[,c("geometry")]
    return(y)
  })
  rm(shapes)
  for(i in 1:4){
    sfs[[i]]$TYPE = cats[i]
    sfs[[i]] = sfs[[i]][,c(2,1)]
  }
  rm(i, cats)
  indicators = rbind(sfs[[1]], sfs[[2]], sfs[[3]], sfs[[4]])
  rm(sfs)
  e = Sys.time()
  print(paste("This took", round(difftime(e, s, units="secs"), 4), "seconds"))
  rm(e, s)
  return(indicators)
}

indicators = read_sf("/Users/aaro9563/Documents/DC_Indicators.shp")
odc = read_sf("/Users/aaro9563/Documents/DC_ODC.shp")
get_em_minutes = function(indicators, odc){
  s = Sys.time()
  minutes = rep(list(list()), 754)
  cats = c("grocery","health","schools","transit")
  for(i in 1:754){
    interest = odc[odc$OriginOID == i,]
    if(nrow(interest) == 0){
      for(j in cats){
        minutes[[i]][[j]] = numeric(0)
      }
    } else{
      types = indicators[interest$Destinat_1,]$TYPE
      for(j in cats){
        minutes[[i]][[j]] = ceiling(interest[which(types == j),]$Total_Time)
      }
    }
  }
  rm(cats, interest, types, i, j)
  e = Sys.time()
  print(paste("This took", round(difftime(e, s, units="secs"), 4), "seconds"))
  rm(e, s)
  return(minutes)
} # Gets minutes from 1 overall call to OD Cost Matrix

# Credit cost for using closest facility to do all routing
approx_creds_all_cf = function(minutes){
  y = lapply(minutes, function(x){
    a = rep(NA, 4)
    for(i in 1:4){
      if(length(x[[i]]) <= 1) a[i]=length(x[[i]])*0.5 else a[i]=1
    }
    return(a)
  })
  z = unlist(y)
  nnr = sum(z == 0)
  nucf1 = sum(z == 0.5)
  nucf2 = sum(z == 1)
  df = data.frame(total_creds = sum(z)+length(z)/4*0.5,
                  sa_creds = length(z)/4*0.5,
                  routing_creds = sum(z),
                  per_use_cf1 = round(nucf1/(length(z)),3),
                  per_use_cf2 = round(nucf2/(length(z)),3),
                  per_no_route = round(nnr/(length(z)),3))
  return(df)
}

# Credit cost for using closest facility when #facilities in SA > 2, find routes otherwise
# Find routes doesn't have cutoff parameter though, so this actually isn't a possibility
approx_creds_part_cf = function(minutes){
  y = lapply(minutes, function(x){
    a = rep(NA, 4)
    for(i in 1:4){
      if(length(x[[i]]) <= 2) a[i]=length(x[[i]])*0.04 else a[i]=1
    }
    return(a)
  })
  z = unlist(y)
  nnr = sum(z == 0)
  znot0 = z[which(z!=0)]
  nufr = sum(znot0-floor(znot0) != 0)
  nucf = sum(znot0-floor(znot0) == 0)
  df = data.frame(total_creds = sum(z)+length(z)/4*0.5,
                  sa_creds = length(z)/4*0.5,
                  routing_creds = sum(z),
                  per_use_cf = round(nucf/(length(z)),3),
                  creds_for_cf = sum(znot0[floor(znot0) != 0]),
                  per_use_fr = round(nufr/(length(z)),3),
                  creds_for_fr = sum(znot0[floor(znot0) == 0]),
                  per_no_route = round(nnr/(length(z)),3))
  return(df)
}

# Credit cost of using ODCM over all sample points and all facilities indicators
approx_creds_overall_odcm = function(minutes, nind){
  df = data.frame(total_creds = length(minutes)*nind*0.0005 + length(minutes)*0.5,
                  sa_creds = length(minutes)*0.5,
                  routing_creds = length(minutes)*nind*0.0005)
}

# Credit cost of using ODCM over each set of facilities indicators and the points for which at least one
# instance of that facilities indicator exists in the service area
approx_creds_segmented_odcm = function(minutes, ng, nh, ns, nt){
  wg = sum(unlist(lapply(minutes, function(x){length(x[[1]]) > 0})))
  wh = sum(unlist(lapply(minutes, function(x){length(x[[2]]) > 0})))
  ws = sum(unlist(lapply(minutes, function(x){length(x[[3]]) > 0})))
  wt = sum(unlist(lapply(minutes, function(x){length(x[[4]]) > 0})))
  df = data.frame(total_creds = 0.0005*(wg*ng + wh*nh + ws*ns + wt*nt) + length(minutes)*0.5,
                  sa_creds = length(minutes)*0.5,
                  routing_creds = 0.0005*(wg*ng + wh*nh + ws*ns + wt*nt),
                  grocery_creds = wg*ng*0.0005,
                  health_creds = wh*nh*0.0005,
                  schools_creds = ws*ns*0.0005,
                  transit_creds = wt*nt*0.0005)
  return(df)
}

rbind(dc = approx_creds_all_cf(dc_minutes),
      chi = approx_creds_all_cf(chi_minutes),
      nyc = approx_creds_all_cf(nyc_minutes),
      la = approx_creds_all_cf(la_minutes)); print("Requires 8 calls to Closest Facility")
rbind(dc = approx_creds_part_cf(dc_minutes),
      chi = approx_creds_part_cf(chi_minutes),
      nyc = approx_creds_part_cf(nyc_minutes),
      la = approx_creds_part_cf(la_minutes)); print("Requires 4 calls to Closest Facility, 4 calls to Find Routes -- THIS DOESN'T WORK")
rbind(dc = approx_creds_overall_odcm(dc_minutes, 267),
      chi = approx_creds_overall_odcm(chi_minutes, 1716),
      nyc = approx_creds_overall_odcm(nyc_minutes, 403),
      la = approx_creds_overall_odcm(la_minutes, 1007)); print("Requires 1 call to OD Cost Matrix")
rbind(dc = approx_creds_segmented_odcm(dc_minutes, 54, 51, 122, 40),
      chi = approx_creds_segmented_odcm(chi_minutes, 506, 42, 1024, 144),
      nyc = approx_creds_segmented_odcm(nyc_minutes, 27, 12, 212, 152),
      la = approx_creds_segmented_odcm(la_minutes, 21, 165, 778, 43)); print("Requires 4 calls to OD Cost Matrix")

