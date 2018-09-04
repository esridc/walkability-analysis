### Step 1: Generating a Sample --------------------

## Induce sample size from city size, street density
get_samp_size = function(city, grid_distance_m = 5*1000*10/60/2, dist_to_road_m = sqrt(2*(5*1000*0.5/60)^2)){
  # grid_distance_m: selected assuming 5km/hr and max-size SA (circular for a 10min walk time)
  #    --> guarantees at least half overlap for each SA
  # dist_to_road_m: selected assuming 5km/hr and min-dist to street associated with is 1 min walk time
  #    --> value is hypotenuse for a right-isosceles with sides lengths = 30 sec walk time
  start_time = Sys.time()
print("Loading libraries...")
  if(!require(sf)){install.packages("sf")}; suppressWarnings(library(sf))
  if(!require(nngeo)){install.packages("nngeo")}; suppressWarnings(library(nngeo))
print("Reading in neighborhoods...")
  neigh = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/Neighborhoods", sep=""), 
                     pattern = ".shp$", full.names = T)
  neighborhoods = read_sf(neigh)
  rm(neigh)
  neighborhoods = st_transform(neighborhoods, 102005)
print("Creating point grid over the city...")
  points = st_make_grid(neighborhoods, cellsize = grid_distance_m, what = "centers")
print("Clipping grid by city...")
  grid_clip = st_intersection(points, neighborhoods)
  grid_clip = st_sf(id = 1:length(grid_clip), geometry = grid_clip)
  rm(points, neighborhoods)
print("Reading in streets...")
  stree = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/City_Streets", sep=""), 
                     pattern = ".shp$", full.names = T)
  streets = read_sf(stree)
  rm(stree)
  streets = st_transform(streets, 102005)
  if(any(is.na(st_dimension(streets)))){
    streets = streets[-which(is.na(st_dimension(streets))),]
  }
print(paste("Finding distance between grid points and streets for", nrow(grid_clip), "grid points..."))
  ndists = c()
  for(i in 1:nrow(grid_clip)){
    nn = st_nn(x=grid_clip[i,], y=streets, k=1, maxdist = dist_to_road_m, returnDist=TRUE, progress=FALSE)
    dmat = nn[[2]]
    rm(nn)
    ndists = append(ndists, dmat[1,1])
    rm(dmat)
    if(i %% 50 == 0){print(paste("Nearest street distance calculated for first", i, "points"))}
  }
  rm(grid_clip, streets, i)
print("Determining sample size...")
  samp_size = sum(!is.na(ndists))
  rm(ndists)
print(paste("Done! Sample size =", samp_size))
  end_time = Sys.time()
  print(paste("Determining sample size took:", end_time - start_time))
  rm(end_time, start_time)
return(samp_size)
}

## Perform stratified regular sampling within a city's neighborhoods based on street coverage
stratified_regular = function(city, samp_size){
  start_time = Sys.time()
print("Loading libraries...")
  if(!require(sf)){install.packages("sf")}; suppressWarnings(library(sf))
  if(!require(sp)){install.packages("sp")}; suppressWarnings(library(sp))
  if(!require(maptools)){install.packages("maptools")}; suppressWarnings(library(maptools))
  if(!require(rgdal)){install.packages("rgdal")}; suppressWarnings(library(rgdal))
print("Reading in neighborhoods...")
  neigh = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/Neighborhoods", sep=""), 
                     pattern = ".shp$", full.names = T)
  neighborhoods = read_sf(neigh)
  rm(neigh)
  neighborhoods = st_transform(neighborhoods, 102005)
  neighborhoods$UNID = 1:nrow(neighborhoods)
print("Reading in streets...")
  stree = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/City_Streets", sep=""), 
                     pattern = ".shp$", full.names = T)
  streets = read_sf(stree)
  rm(stree)
  streets = st_transform(streets, 102005)
  if(any(is.na(st_dimension(streets)))){
    streets = streets[-which(is.na(st_dimension(streets))),]
  }
print("Intersecting streets and neighborhoods...")
  int = suppressWarnings(st_intersection(streets, neighborhoods))
print("Obtaining number of samples in each neighborhood by percent coverage of number of streets...")
  street_num = c()
  for(i in 1:nrow(neighborhoods)){
    street_num = append(street_num, sum(int$UNID == i))
  }
  percent_streets = street_num / nrow(int)
  rm(i, street_num)
  samp_by_neighborhood = ceiling(percent_streets * samp_size)
  rm(percent_streets)
  sbn_df = data.frame("neighb_ID" = neighborhoods$UNID, "samp_size" = samp_by_neighborhood)
print("Sample sizes - overall and per-neighborhood - calculated")
print(paste("Expected sample size:", samp_size))
print(paste("Actual sample size:  ", sum(sbn_df[,2])))
print("Per neighborhood sample sizes:")
print(sbn_df)
print("Converting to requisite form for sampling...")
  neighborhoods = as(neighborhoods, "Spatial")
  streets = as(streets, "Spatial")
print(paste("Obtaining the sample points and snapping to streets for", nrow(neighborhoods), "neighborhoods..."))
  raw = spsample(neighborhoods[neighborhoods$UNID == 1,], n = samp_by_neighborhood[1], type = "regular")
  samp = snapPointsToLines(raw, streets)
  for(i in 2:nrow(neighborhoods)){
    raw = spsample(neighborhoods[neighborhoods$UNID == i,], n = samp_by_neighborhood[i], type = "regular")
    raw = snapPointsToLines(raw, streets)
    samp = rbind(samp, raw)
    if(i %% 10 == 0){print(paste("Samples obtained for the first", i, "neighborhoods"))}
  }
  rm(i, raw, samp_by_neighborhood)
print("Writing sample shapefile...")
  sampdir = "/Users/aaro9563/Documents/Walk_Final/"
  if(!file.exists(sampdir)){dir.create(sampdir)}
  suppressWarnings(writeOGR(obj = samp, dsn = sampdir, layer = paste(city, "_Stratified_Regular", sep=""), driver = "ESRI Shapefile",
                            overwrite_layer = TRUE))
  rm(samp, sampdir)
print(paste("Done! Refer to the ", city, "_Stratified_Regular shapefile for service area creation", sep=""))
  end_time = Sys.time()
  print(paste("Sampling itself took:", end_time - start_time))
  rm(end_time, start_time)
}

## Perform the whole sampling algorithm
walkability_sampling = function(city){
  s = Sys.time()
  print("INDUCING SAMPLE SIZE FOR SPECIFIED CITY")
  ss = get_samp_size(city=city)
  print("PERFORMING STRATIFIED-REGULAR SAMPLING (BY NEIGHBORHOOD)")
  stratified_regular(city=city, samp_size = ss)
  e = Sys.time()
  print(paste("The entire sampling process took:", e-s))
  rm(s, e, ss)
}

### Step 2: Generating Service Areas (Go to Arc to do this!) --------------------

# Once service areas are created, save the shapefile to: Documents/Walk_Final/AOI_Service_Areas.shp
# Could use the function written in "ArcGIS API for Raster.py"

# If the sample size is > 1000, here's a format to bind batch-process-generated service areas together
# I.e. this is a general framework for how you can write a series of 1000 point service areas (created in Arc) into one combined file (like we need)
# In this example, la01 is the servicea areas for points 1-1000, la12 is the service areas for points 1001-2000, ...

la01 = read_sf("/Users/aaro9563/Documents/LA_Service_Areas_01.shp")
la12 = read_sf("/Users/aaro9563/Documents/LA_Service_Areas_12.shp")
la23 = read_sf("/Users/aaro9563/Documents/LA_Service_Areas_23.shp")
la34 = read_sf("/Users/aaro9563/Documents/LA_Service_Areas_34.shp")
la45 = read_sf("/Users/aaro9563/Documents/LA_Service_Areas_45.shp")

la01 = la01[order(la01$FacilityOI),]
la12 = la12[order(la12$FacilityOI),]
la23 = la23[order(la23$FacilityOI),]
la34 = la34[order(la34$FacilityOI),]
la45 = la45[order(la45$FacilityOI),]

la12$FacilityOI = la12$FacilityOI + 1
la23$FacilityOI = la23$FacilityOI + 1
la34$FacilityOI = la34$FacilityOI + 1
la45$FacilityOI = la45$FacilityOI + 1

la = rbind(la01, la12, la23, la34, la45)
st_write(la, "/Users/aaro9563/Documents/Walk_Final/LA_Service_Areas.shp")

rm(la, la01, la12, la23, la34, la45)

### Step 3: Routing to Facilities Indicators --------------------

# Making a graph... help from Barry Rowlingson RPubs doc on route finding
# Thanks Barry! http://rpubs.com/geospacedman/routing
create_graph_from_lines = function(sp_lines) {
  if(!require(sp)){install.packages("sp")}; suppressWarnings(library(sp))
  if(!require(igraph)){install.packages("igraph")}; suppressWarnings(library(igraph))
  if(!require(rgeos)){install.packages("rgeos")}; suppressWarnings(library(rgeos))
  # Intersect lines with each other to obtain all the verticies
  inters = gIntersection(sp_lines, sp_lines)
  # Get coordinates of the segments resulting from the intersection
  edge = do.call(rbind, lapply(inters@lines[[1]]@Lines, function(x){
    as.vector(t(x@coords))
    })
  )
  rm(inters)
  # For edge weights, use euclidean distance as an approximation of actual distance
  # Should work well considering we're in cities and lines will be short (regardless of curviness)
  lens = sqrt((edge[,1] - edge[,3])^2 + (edge[,2] - edge[,4])^2) 
  # Organize the coordinates on which to build the graph
  # Leave as character vectors so they aren't interpreted as vertex IDs
  # Instead, they'll be names (which we'll use later)
  from = paste(edge[,1], edge[,2], sep=",")
  to = paste(edge[,3], edge[,4], sep=",")
  build = cbind(from, to)
  rm(edge, from, to)
  # Make the graph, assign the edge weights as the lengths (see above)
  citygraph = graph.edgelist(build, directed = FALSE)
  E(citygraph)$weight = lens
  rm(lens, build)
  # Extract the coordinates from the names of the graph elements, and use to create the vertices
  verts = do.call(rbind, strsplit(V(citygraph)$name, ","))
  V(citygraph)$x = as.numeric(verts[,1])
  V(citygraph)$y = as.numeric(verts[,2])
  rm(verts)
  # We out
  return(citygraph)
}

# Extract shortest path distance along streets from sample point to indicator
# Thanks again Barry for the framework^^
shortest_route_along_lines = function(city_graph, sample_point, indicator_point){
  # Get nearest vertices to point, and distances to those vertices
  # Again, being in cities, those distances ~should~ be small
  # Here, we'll add the euclidean distances as "error" between actual point and vertex
  if(!require(igraph)){install.packages("igraph")}; suppressWarnings(library(igraph))
  if(!require(FNN)){install.packages("FNN")}; suppressWarnings(library(FNN))
  verts = cbind(V(city_graph)$x, V(city_graph)$y)
  from_nn = get.knnx(verts, sample_point, 1)
  to_nn = get.knnx(verts, indicator_point, 1)
  rm(verts)
  from_vert = from_nn$nn.index[1,1]
  from_error = from_nn$nn.dist[1,1]
  rm(from_nn)
  to_vert = to_nn$nn.index[1,1]
  to_error = to_nn$nn.dist[1,1]
  rm(to_nn)
  # Get the shortest path distance
  p = shortest.paths(city_graph, from_vert, to_vert)
  path_dist = p[1,1]
  rm(p, from_vert, to_vert)
  est_dist = path_dist + from_error + to_error
  rm(path_dist, from_error, to_error)
  return(est_dist)
}

# Extract them minutes
minute_extraction_by_graph = function(city){
  s = Sys.time()
print("EXTRACTING WALK TIME DISTANCES FROM SAMPLES TO WALKABLE INDICATORS")
print("Loading libraries...")
  if(!require(sf)){install.packages("sf")}; suppressWarnings(library(sf))
  if(!require(sp)){install.packages("sp")}; suppressWarnings(library(sp))
  if(!require(rgdal)){install.packages("rgdal")}; suppressWarnings(library(sf))
print("Reading in data...")
  samplespath = paste("/Users/aaro9563/Documents/Walk_Final/", city, "_Stratified_Regular.shp", sep="")
  servarspath = paste("/Users/aaro9563/Documents/Walk_Final/", city, "_Service_Areas.shp", sep="")
  streetspath = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/City_Streets", sep=""), pattern = ".shp$", full.names=TRUE)[1]
  grocerypath = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/DUC Input/Grocery", sep=""), pattern = ".shp$", full.names=TRUE)[1]
  healthspath = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/DUC Input/Health", sep=""), pattern = ".shp$", full.names=TRUE)[1]
  schoolspath = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/DUC Input/Schools", sep=""), pattern = ".shp$", full.names=TRUE)[1]
  transitpath = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/DUC Input/Transit", sep=""), pattern = ".shp$", full.names=TRUE)[1]
  samplepoints = readOGR(samplespath, verbose=FALSE)
  serviceareas = read_sf(servarspath)
  citysstreets = read_sf(streetspath)
  if(any(is.na(st_dimension(citysstreets)))){
    citysstreets = citysstreets[-which(is.na(st_dimension(citysstreets))),]
  }
  citysstreets = as(citysstreets, "Spatial")
  grocerystore = read_sf(grocerypath)
  if(any(is.na(st_dimension(grocerystore)))){
    grocerystore = grocerystore[-which(is.na(st_dimension(grocerystore))),]
  }
  healthcenter = read_sf(healthspath)
  if(any(is.na(st_dimension(healthcenter)))){
    healthcenter = healthcenter[-which(is.na(st_dimension(healthcenter))),]
  }
  publicschool = read_sf(schoolspath)
  if(any(is.na(st_dimension(publicschool)))){
    publicschool = publicschool[-which(is.na(st_dimension(publicschool))),]
  }
  transitstops = read_sf(transitpath)
  if(any(is.na(st_dimension(transitstops)))){
    transitstops = transitstops[-which(is.na(st_dimension(transitstops))),]
  }
  rm(streetspath, samplespath, servarspath, grocerypath, healthspath, schoolspath, transitpath)
print("Projecting data to USA Equidistant Conic...")
  usaeqc = CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
  citysstreets = spTransform(citysstreets, usaeqc)
  samplepoints = spTransform(samplepoints, usaeqc)
  serviceareas = st_transform(serviceareas, 102005)
  grocerystore = st_transform(grocerystore, 102005)
  healthcenter = st_transform(healthcenter, 102005)
  publicschool = st_transform(publicschool, 102005)
  transitstops = st_transform(transitstops, 102005)
  rm(usaeqc)
print("Creating graph from city streets...")
  citygraph = create_graph_from_lines(sp_lines = citysstreets)
  rm(citysstreets)
print("Organizing data...")
  indicators = list("grocery"=grocerystore,
                    "health"=healthcenter,
                    "schools"=publicschool,
                    "transit"=transitstops)
  rm(grocerystore, healthcenter, publicschool, transitstops)
  n = nrow(samplepoints)
  minutes = rep(list(list(c(),c(),c(),c())), n)
print(paste("Extracting walk times to walkable indicators for", n, "sample points..."))
  for(i in 1:n){
    pointcoords = data.frame("x"=unname(samplepoints@coords[i,1]),
                             "y"=unname(samplepoints@coords[i,2]))
    sa = serviceareas[serviceareas$FacilityOI==i,]
    for(j in 1:4){
      ind = indicators[[j]]
      inters = suppressWarnings(st_intersection(ind, sa)) #all it tells us is we want spatially consistent attributes, which for us doesn't matter
      rm(ind)
      if(nrow(inters) == 0){
        rm(inters)
        minutes[[i]][[j]] = numeric(0)
      } else{
        inters = as(inters, "Spatial")
        indiccoords = data.frame("x"=unname(inters@coords[,1]),
                                 "y"=unname(inters@coords[,2]))
        rm(inters)
        vec = c()
        for(k in 1:nrow(indiccoords)){
          d = shortest_route_along_lines(city_graph = citygraph,
                                         sample_point = pointcoords,
                                         indicator_point = indiccoords[k,])
          vec = append(vec,d)
        }
        rm(indiccoords, d, k)
        vec = vec/1000*60/5
        vec = ceiling(vec)
        vec = vec[vec <= 10]
        minutes[[i]][[j]] = vec
        rm(vec)
      }
      rm(j)
    }
    if(i %% 20 == 0){print(paste("Walk-time distances obtained for the first", i, "sample points"))}
    rm(pointcoords, sa, i)
  }
print("Done! Minutes successfully extracted to a list (indexed by sample point ID) of lists (indexed as grocery, health, school, transit).")
  rm(citygraph)  
  e = Sys.time()
  print(paste("This process took:", e-s))
  rm(e, s)
  return(minutes)
}

### Step 4: Obtaining Non-Accessibility Qualities --------------------

# Calculate the percentages of area indicators in each service area
swing_input = function(city){
  s = Sys.time()
print("EXTRACTING CHARACTERISTICS OF SAMPLE POINTS' WALKABLE AREA")
print("Loading libraries...")
  if(!require(sf)){install.packages("sf")}; suppressWarnings(library(sf))
print("Reading in the service area polygons...")
  service_area = read_sf(paste("/Users/aaro9563/Documents/Walk_Final/", city, "_Service_Areas.shp", sep=""))
  service_area = st_transform(service_area, 102005)
  service_area = service_area[order(service_area$FacilityOI),]
print("Reading in city characteristic files...")
  direcs = list.dirs(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/Decay Parameter Input", sep=""))[-1]
  shp = c()
  for(folder in direcs){
    shp = append(shp, list.files(folder, pattern = ".shp$", full.names = T))
  }
  per_table = data.frame("point"=as.numeric(service_area$FacilityOI))
  ord = c("first", "second", "third", "fourth")
  rm(folder, direcs)
print(paste("Calculating percentages of each city characteristic in each of", nrow(service_area), "service areas..."))
  for(doc in shp){
    myshp = read_sf(doc)
    myshp = st_transform(myshp, 102005)
    n = nrow(myshp)
    if(any(is.na(st_dimension(myshp)))){
      myshp = myshp[-which(is.na(st_dimension(myshp))),]
    }
    counts = lengths(st_covers(service_area, myshp))
    counts = unname(counts)
    per_table = cbind(per_table, counts/n)
    rm(n, myshp, counts)
    print(paste("Percentages calculated for the", ord[which(shp == doc)], "of four city characteristic files"))
  }
  rm(doc, shp, ord, service_area)
  names(per_table) = c("point", "crashes", "crime", "historic", "trees")
print("Done! Per-service area percentages successfully extracted to a data frame (rows: service areas; columns: city characteristics")
  e = Sys.time()
  print(paste("This process took:", e-s))
  rm(e,s)
  return(per_table)
}

### Step 5: Applying the Walkability model --------------------

# Minimax normalizer -- keeps things on a 0-1 scale
minimax = function(x){
  (x-min(x))/(max(x)-min(x))
}

# Solving for steepness [for proximity rank weighting function]
# (numerical solve, given a desired weight for the (n_c)th or (n_c + 1)th point in DUC c in a service area 
# [where n_c is the user selected desired number of points of DUC c in a service area]
# use Python Bridge because R solvers don't get us close enough to 0)
k_est = function(fl = 100, dn = 1, bw = 0.2, ...){
  python.assign("l", fl)
  python.assign("n_c", dn)
  python.assign("p", bw)
  python.exec('
import numpy
from scipy.optimize import fsolve
if p < 0.5:
  func = lambda k : (1/(1+numpy.exp(k*(n_c + 1 - n_c - 0.5))) - 1/(1+numpy.exp(k*(l - n_c - 0.5)))) / (p*(1/(1+numpy.exp(k*(1 - n_c - 0.5))) - 1/(1+numpy.exp(k*(l - n_c - 0.5))))) - 1
  k_init = 2*numpy.log((1 - p)/p)
  k_sol = fsolve(func, k_init)
  k_sol = k_sol[0]
elif p > 0.5 and n_c != 1:
  func = lambda k : (1/(1+numpy.exp(k*(n_c - n_c - 0.5))) - 1/(1+numpy.exp(k*(l - n_c - 0.5)))) / (p*(1/(1+numpy.exp(k*(1 - n_c - 0.5))) - 1/(1+numpy.exp(k*(l - n_c - 0.5))))) - 1
  k_init = 2*numpy.log((1 - p)/p)
  k_sol = fsolve(func, k_init)
  k_sol = -k_sol[0]
elif p > 0.5 and n_c == 1:
  p = 1-p
  func = lambda k: (1 / (1 + numpy.exp(k * (n_c + 1 - n_c - 0.5))) - 1 / (1 + numpy.exp(k * (l - n_c - 0.5)))) / (p * (1 / (1 + numpy.exp(k * (1 - n_c - 0.5))) - 1 / (1 + numpy.exp(k * (l - n_c - 0.5))))) - 1
  k_init = 2*numpy.log((1 - p)/p)
  k_sol = fsolve(func, k_init)
  k_sol = k_sol[0]
')
  k = python.get('k_sol')
  # print(paste("Steepness parameter =", k))
  return(k)
}

# Proximity rank weighting function -- weights based on distance rank from an origin point
close = function(fun_len = 100, desired_number = 1, boundary_weight = .2, steepness = NULL, ...){
  l = seq(1, fun_len, by=1)
  if(is.null(steepness)){
    k = k_est(fl = fun_len, dn = desired_number, bw = boundary_weight)
    steep = k
  } else {
    steep = steepness
    print(paste("Steepness parameter =", steep))
  }
  close_df = data.frame(l, w = 1/(1+exp(steep*(l - desired_number - 0.5))))
  close_df$w = minimax(close_df$w)
  normalizer = sum(close_df$w[1:desired_number])
  # p = ggplot(data = close_df) +
  #   geom_line(mapping = aes(x = l, y = w)) +
  #   xlim(c(1, desired_number + 10)) +
  #   labs(x = "Point Rank (by Time to Origin)",
  #        y = "Weight",
  #        title = "Decay Function for Closeness Weighting",
  #        subtitle = paste("Desired Number Walkable = ", desired_number, ", Steepness = ", round(steep,4), sep = ""))
  # print(p)
  return(list(close_df = close_df, normalizer = normalizer))
}

# Area characteristics parameter -- +- for a score based on area characteristics
area_swing = function(crashes = 0, crime = 0, historic = 0, trees = 0, ...){
  swing = 0 + 100*(- crashes - crime + historic + trees)
  # print(paste("Area Swing Parameter =", round(swing,4)))
  return(swing)
}

# Distance weighting function -- weights based on walk time distance from an origin point
dis = function(walkable_distance = 10, ...){
  d = seq(1, walkable_distance, by = 1)
  dis_df = data.frame(d, w = 1-d)
  dis_df$w = minimax(dis_df$w)*0.5 + 0.5
  # p = ggplot(data = dis_df) +
  #   geom_line(mapping = aes(x = d, y = w)) +
  #   labs(x = "Distance (in Minutes)",
  #        y = "Weight",
  #        title = "Decay Function for Distance Weighting",
  #        subtitle = paste("Walkable Distance = ", walkable_distance, " minutes", sep = ""))
  # print(p)
  return(dis_df)
}

# Category walk score function -- walk score for a particular facilities indicator
weighter = function(minutes, ...){
  if(length(minutes) == 0){
    walk = 0
    # print(paste("Category Walk Score =", walk, "--- no points in a walkable distance of the origin"))
  } else {
    in_area = data.frame(count = 1:length(minutes), min = sort(minutes))
    dis_decay = dis(...)
    close_decay = close(...)
    dweight = c()
    pweight = c()
    for(i in in_area$min){
      dweight = append(dweight, dis_decay[dis_decay$d == i, 2])
    }
    for(i in in_area$count){
      pweight = append(pweight, close_decay$close_df[close_decay$close_df$l == i, 2])
    }
    weight = cbind(in_area, dweight, pweight)
    names(weight) = c("count","distance","dist_weight","count_weight")
    walk_first = sum(weight$dist_weight * weight$count_weight)
    if(walk_first/close_decay$normalizer > 1){
      walk = 1
    } else {
      walk = walk_first/close_decay$normalizer
    }
    # print("Table of Points, Distances, and Weights:")
    # print(weight)
    # print(paste("Category Walk Score = ", walk, " (from unbounded sum ", walk_first, " divide best-possible-scenario ", close_decay$normalizer, ")", sep = ""))
  }
  return(walk)
}

# Get a vector of walkability scores
walkability = function(minutes_list, swing_inputs, ...){
  s = Sys.time()
  print("CALCULATING WALKABILITY SCORES FOR SAMPLE POINTS")
  print("Loading libraries...")
  if(!require(rPython)){install.packages("rPython")}; suppressWarnings(library(rPython))
  print("Defining inputs...")
  inputs = data.frame(category_weight   = c(.25, .25, .25, .25),
                      walkable_distance = c(10 , 10 , 10 , 10),
                      desired_number    = c(1  , 1  , 1  , 1))
  rownames(inputs) = c("grocery", "health", "schools", "transit")
  bound = 0.2
  s_df = swing_inputs   # swing_inputs must come from the saved output of the "swing_input" function
  print(paste("Calculating walk scores for", nrow(s_df), "sample points..."))
  walk_scores = rep(NA, length(minutes_list))
  for(i in 1:length(walk_scores)){
    cs = rep(NA, nrow(inputs))
    for(j in 1:nrow(inputs)){
      w = inputs$category_weight[j] * 100
      cs[j] = w * weighter(minutes = minutes_list[[i]][[j]], 
                           walkable_distance = inputs$walkable_distance[j], 
                           desired_number = inputs$desired_number[j], boundary_weight = bound)
    }
    swing = area_swing(crashes = s_df[i,"crashes"], crime = s_df[i,"crime"], historic = s_df[i,"historic"], trees = s_df[i,"trees"])
    indscore = sum(cs)
    if(abs(swing) > indscore){
      adder = indscore*sign(swing)
    } else{
      adder = swing
    }
    total_walk = indscore + adder
    if(total_walk < 0){
      total_walk = 0
    }
    walk_scores[i] = total_walk
  }
  print("Done! Walk scores successfully extracted to a vector (indexed by sample point ID)")
  e = Sys.time()
  print(paste("This process took:", e-s))
  rm(e,s)
  return(walk_scores)
}


### Step 6: (Write the walk scores to the sample points for) Generating the Raster --------------------

# Write that file
walk_enriching = function(city, walkscore_vector){
  s = Sys.time()
print("WRITING WALK SCORES TO SAMPLE POINTS")
print("Loading libraries...")
  if(!require(rgdal)){install.packages("rgdal")}; suppressWarnings(library(rgdal))
print("Writing shapefile...")
  samp = readOGR(paste("/Users/aaro9563/Documents/Walk_Final/", city, "_Stratified_Regular.shp", sep=""),
                 verbose = FALSE)
  samp@data$walk = walkscore_vector
  writeOGR(obj = samp, 
           dsn = "/Users/aaro9563/Documents/Walk_Final",
           layer = paste(city, "_Walk_Enriched", sep=""),
           driver = "ESRI Shapefile",
           overwrite_layer = TRUE)
print("Done! Shapefile of sample points with walkability scores attribute successfully written.")
  e = Sys.time()
  print(paste("This process took:", e-s))
  rm(e,s)
}



### THIS FUNCTION RUNS THE ALGORITHM -- STEPS 1-6 --------------------

# This function runs the whole algorithm. These are the inputs:
# AOI: The city. Needs to be exactly the same as whatever you're calling the folder where your data is
#   ex. my Manhattan data was in a folder called NYC. So AOI="NYC"
# got_sample: Have you already created the sample? If yes, change this to anything but "no", and sample generation will be skipped
# got_sa: Have you already created service areas? If yes, change this to anything but "no", and you won't be prompted to go create SA
# minutes: if you've already created the list of list of walk-time distances, enter it here
# swingtable: if you've already created the dataframe of area indicator percentages, enter it here,
# walkscores: if you've already created the vector of walkability scores, enter it here
# write: Have you already written the walkability scores to the sample? If yes, change this to anything but "yes"
#   ex. Note that if you're changing this parameter, you should probably just run whatever you want to run individually

# IF YOU'RE USING NEW DATA, IT'S IMPORTANT THAT IT FOLLOWS THE SAME FILE STRUCTURE THAT I'VE SET UP ALREADY
# IF IT DOESN'T, THIS WILL NOT RUN
# (well, I guess it would if you went an changed all the read and write locations...)
# Refer to my files to see what that looks like. Something like this:
# (Note: I'm only including what you absolutely need, my files include some other things that aren't necessary for this)
# > Documents
# >>> 01 Walkability
# >>>>>> Some city
# >>>>>>>>> City General
# >>>>>>>>>>>> City_Streets -- this folder should include only the streets shapefile
# >>>>>>>>>>>> Neighborhoods  -- this folder should include only the neighborhoods shapefile
# >>>>>>>>> DUC Input
# >>>>>>>>>>>> Grocery-- this folder should include only the grocery shapefile
# >>>>>>>>>>>> Health -- this folder should include only the health shapefile
# >>>>>>>>>>>> Schools -- this folder should include only the schools shapefile
# >>>>>>>>>>>> Transit -- this folder should include only the transit shapefile
# >>>>>>>>> Decay Parameter Input
# >>>>>>>>>>>> Crashes -- this folder should include only the crashes shapefile
# >>>>>>>>>>>> Crime -- this folder should include only the crime shapefile
# >>>>>>>>>>>> Historic Sites -- this folder should include only the historic sites shapefile
# >>>>>>>>>>>> Street Trees -- this folder should include only the street trees shapefile

# Okay, here's the function. 
# Again, this will complete the process end to end, with the exception of service area creation
# Please refer to step 2 to see how/where service areas should be saved!
# This function will return the minutes list, area-indicators table, and walk score vector
# I'd recommend saving newly produced versions of these as .Rdata files somewhere so you don't have to keep re-creating them

create_walkability_raster = function(AOI, got_sample="no", got_sa="no", 
                                     minutes, swingtable, walkscores, write="yes"){
  s = Sys.time()
print(paste("WELCOME! WE'RE GOING TO CREATE THE WALKABILITY SCORE FOR", AOI))
  if(got_sample == "no"){
    walkability_sampling(city=AOI)
  } else{
    print("Sample already created")
  }
  if(got_sa == "no"){
    print("GO TO ARC AND CREATE SERVICE AREAS")
    print("Save as '/Users/aaro9563/Documents/Walk_Final/AOI_Service_Areas.shp'")
    print("Once you've done this, return here and type 'yes' then press 'enter'")
    sa_done = readline(prompt = "Have service areas been created and saved? ")
  } else{
    print("Service areas already created")
    sa_done = 'yes'
  }
  if(sa_done == 'yes'){
    if(missing(minutes)){
      minutes = minute_extraction_by_graph(city=AOI)
    } else{
      print("Minutes list for sample given")
    }
    if(missing(swingtable)){
      swingtable = swing_input(city=AOI)
    } else{
      print("Swing inputs data frame for sample given")
    }
    if(missing(walkscores)){
      walkscores = walkability(minutes_list=minutes, swing_inputs=swingtable)
    } else{
      print("Walk score vector for sample given")
    }
    if(write == "yes"){
      walk_enriching(city=AOI, walkscore_vector=walkscores)
      print(paste("PROCESS COMPLETE -- SAMPLE AND ASSOCIATED WALKSCORES GENERATED FOR", AOI))
      print("Walk enriched shapefile written to '/Users/aaro9563/Documents/Walk_Final/AOI_Walk_Enriched.shp'")
      print("Proceed to ArcGIS to create the raster, using IDW with parameters search=2, variable=4")
    } else{
      print("Shapefile writing not requested -- output is walk score vector")
      print(paste("PROCESS COMPLETE -- SAMPLE AND ASSOCIATED WALKSCORES GENERATED FOR", AOI))
    }
  }
  e = Sys.time()
  print(paste("The entire raster process, including time spent to create service areas, took:", e-s))
  rm(e,s)
  return(list(minutes, swingtable, walkscores))
}




### Obtaining summary stats for a city APPROXIMATION --------------------

# DISCLAIMER: THIS IS JUST AN APPROXIMATION FUNCTION. IT IS NOT THE REAL THING
# In real time, we'll be querying the raster, and using percentiles
# This gives point score (of closest sample point in neighborhood), neighborhood median, and city median
# Basically, it geocodes a point, finds its neighborhood, and then gives all of the above
# Using geocoder="google" limits you to 2500 queries/day... prolly not an issue, but if it is, use geocoder="dsk" instead
# But, the "dsk" geocoder is not nearly as accurate
# Note that sometimes google geocoder can be finicky. If it doesn't go, just try running it again.

# Get point, neighborhood, city scores given an address
summary_scores = function(address, city, geocoder="google"){
  s = Sys.time()
print("OBTAINING SUMMARY SCORE INFORMATION FOR THE PROVIDED ADDRESS")
print("Loading libraries...")
print("Loading libraries...")
  if(!require(sf)){install.packages("sf")}; suppressWarnings(library(sf)); print("sf loaded")
  if(!require(ggmap)){install.packages("ggmap")}; suppressWarnings(library(ggmap)); print("ggmap loaded")
  if(!require(nngeo)){install.packages("nngeo")}; suppressWarnings(library(nngeo)); print("nngeo loaded")
print("Geocoding address...")
  gcadd = geocode(address, output = "latlon", source = geocoder)
  gcadd$ID = 1
  gcadd = st_as_sf(gcadd, coords = c("lon", "lat"), crs = 4326)
  gcadd = st_transform(gcadd, 102005)
print("Finding neighborhood of address...")
  neigh = list.files(paste("/Users/aaro9563/Documents/01 Walkability/", city, "/City General/Neighborhoods", sep=""), 
                     pattern = ".shp$", full.names = T)
  neighborhood = read_sf(neigh)
  rm(neigh)
  neighborhood = st_transform(neighborhood, 102005)
  yn = unlist(lapply(st_contains(neighborhood, gcadd), function(x){length(x)}))
  neighborhood = neighborhood[which(yn == 1),]
  rm(yn)
print("Finding all samples in that neighborhood...")
  sample_points = read_sf(paste("/Users/aaro9563/Documents/Walk_Final/", city, "_Walk_Enriched.shp", sep=""))
  sample_points = st_transform(sample_points, 102005)
  consideration_points = suppressWarnings(st_intersection(sample_points, neighborhood))
  rm(neighborhood)
print("1. Getting est. address score (score of nearest sample point)...")
  nearest_id = st_nn(gcadd, consideration_points, k=1)[[1]][1]
  point_score = consideration_points[nearest_id,]$walk
  rm(nearest_id)
print("2. Getting est. neighborhood score (median score of sample points in neighborhood)...")
  neighborhood_median = median(consideration_points$walk)
  rm(consideration_points)
print("3. Getting est. city score (median score of all sample points in city)...")
  city_median = median(sample_points$walk)
  scores_df = data.frame("point"=point_score, "neighborhood"=neighborhood_median, "city"=city_median)
  rm(point_score, neighborhood_median, city_median, sample_points)
print("Done! All summary scores successfully saved")
  e = Sys.time()
  print(paste("This process took:", e-s))
  rm(e,s)
  return(scores_df)
}

# Example use of this using some points in DC
whitehouse = summary_scores("1600 Pennsylvania Ave NW, Washington DC", "DC")
cathedral = summary_scores("3101 Wisconsin Ave NW, Washington DC", "DC")
zoo = summary_scores("3001 Connecticut Ave NW, Washington DC", "DC")
art = summary_scores("600 Constitution Ave NW, Washington DC", "DC")
capitol = summary_scores("101 E Capitol St NE, Washington DC", "DC")
natspark = summary_scores("1500 S Capitol St SE, Washington DC", "DC")
df = rbind("White House"=whitehouse,
           "US Capitol"=capitol,
           "Natl. Gallery of Art"=art,
           "Natl. Cathedral"=cathedral,
           "Natl. Zoo"=zoo,
           "Nationals Park"=natspark)
df = round(df, 2)



