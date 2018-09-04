read_decay_fit_data = function(){
  s = Sys.time()
print("Loading libraries...")
  if(!require(data.table)){install.packages("data.table")}; library(data.table)
  if(!require(bit64)){install.packages("bit64")}; library(bit64)
print("Reading in full data...")
  nhts = fread("/Users/aaro9563/Documents/O1 Walkability/NHTS Data/Ascii/DAYV2PUB.CSV")
  dat = list()
print("Subsetting data for all non-transit indicators...")
  nhts_not_transit = nhts[nhts$TRPTRANS == 23 & nhts$TRVLCMIN != -9 & nhts$USEPUBTR != 01,] #walking trips
  nhts_not_transit = nhts_not_transit[,c("TRVLCMIN", "WHYTO", "URBANSIZE", "HH_CBSA", "HHSTATE")]
  dat$all_but_transit = nhts_not_transit
  rm(nhts_not_transit)
print("Defining trip purpose codes for all non-transit indicators...")
  use_codes = data.frame(rbind(
    c(21, "school"),
    c(30, "health"),
    c(41, "grocery")
  ))
  names(use_codes) = c("Code", "Purpose")
  dat$use_codes = use_codes
  rm(use_codes)
print("Subsetting data for transit indicator...")
  nhts_transit = nhts[nhts$TRACC1 == 23 & nhts$PUBTRANS == 01,]
  nhts_transit = nhts_transit[,c("TRACCTM", "URBANSIZE", "HH_CBSA", "HHSTATE")]
  dat$just_transit = nhts_transit
  rm(nhts_transit, nhts)
print("Defining city metro area codes...")
  city_codes = data.frame(rbind(
    c(12060, "Atlanta, GA"),
    c(12420, "Austin, TX"),
    c(12580, "Baltimore, MD"),
    c(13820, "Birmingham, AL"),
    c(14460, "Boston, MA"),
    c(15380, "Buffalo, NY"),
    c(16740, "Charlotte, NC"),
    c(16980, "Chicago, IL"),
    c(17140, "Cincinnati, OH"),
    c(17460, "Cleveland, OH"),
    c(18140, "Columbus, OH"),
    c(19100, "Dallas-Fort Worth, TX"),
    c(19740, "Denver, CO"),
    c(19820, "Detroit, MI"),
    c(25540, "Hartford, CT"),
    c(26420, "Houston, TX"),
    c(26900, "Indianapolis, IN"),
    c(27260, "Jacksonville, FL"),
    c(28140, "Kansas City, MO"),
    c(29820, "Las Vegas, NV"),
    c(31100, "Los Angeles, CA"),
    c(31140, "Louisville, KY"),
    c(32820, "Memphis, TN"),
    c(33100, "Miami, FL"),
    c(33340, "Milwaukee, WI"),
    c(33460, "Minneapolis-St. Paul, MN"),
    c(34980, "Nashville, TN"),
    c(35380, "New Orleans, LA"),
    c(35620, "New York, NY"),
    c(36420, "Oklahoma City, OK"),
    c(36740, "Orlando, FL"),
    c(37980, "Philadelphia, PA"),
    c(38060, "Phoenix, AZ"),
    c(38300, "Pittsburgh, PA"),
    c(38900, "Portland, OR"),
    c(39300, "Providence, RI"),
    c(39580, "Raleigh, NC"),
    c(40060, "Richmond, VA"),
    c(40140, "Riverside, CA"),
    c(40380, "Rochester, NY"),
    c(40900, "Sacramento, CA"),
    c(41180, "St. Louis, MO"),
    c(41620, "Salt Lake City, UT"),
    c(41700, "San Antonio, TX"),
    c(41740, "San Diego, CA"),
    c(41860, "San Francisco, CA"),
    c(41940, "San Jose, CA"),
    c(42660, "Seattle, WA"),
    c(45300, "Tampa, FL"),
    c(47260, "Virginia Beach, VA"),
    c(47900, "Washington, DC")
  ))
  names(city_codes) = c("Code", "Metro")
  dat$city_codes = city_codes
  rm(city_codes)
print("Done! Data successfully compiled.")
  e = Sys.time()
  difftime(e, s, units="mins")
  rm(e, s)
  return(dat)
}

distance_weights_by_city = function(metro_area, nhts_list_of_tables, ss){
  s = Sys.time()
print("Extracting AOI information from parameter...")
  cc = nhts_list_of_tables$city_codes
  city = cc[cc$Metro == metro_area, "Code"]
  state = strsplit(metro_area, ", ")[[1]][2]
  rm(cc)
print("Constructing weights for walk-time distances to all non transit indicators...")
print("(based on the AOI, where granularity is a function of NHTS samples occuring in the AOI)")
  rest = nhts_list_of_tables$all_but_transit
  n_city = sum(rest$HH_CBSA == city)
  n_state = sum(rest$HHSTATE == state)
  weights_list = list()
  for(i in c("grocery", "health", "school")){
    purpose = nhts_list_of_tables$use_codes
    purpose = purpose[purpose$Purpose == i, "Code"]
    if(n_city >= ss){
      print(paste("Weights for", i, "will be based on trips in the provided metro area"))
      times = rest[rest$HH_CBSA == city & rest$WHYTO == purpose, "TRVLCMIN"]
      ecdf_of_times = ecdf(times)
      correction = ecdf_of_times(1)
      weights_list$i = list(ecdf_of_times, correction)
      rm(times, ecdf_of_times, correction)
    } else{
      if(n_state >= ss){
        print(paste("Weights for", i, "will be based on trips in the provided metro area's state"))
        times = rest[rest$HHSTATE == state & rest$WHYTO == purpose, "TRVLCMIN"]
        ecdf_of_times = ecdf(times)
        correction = ecdf_of_times(1)
        weights_list$i = list(ecdf_of_times, correction)
        rm(times, ecdf_of_times, correction)
      } else{
        print(paste("Weights for", i, "will be based on trips in all available locations"))
        times = rest[rest$WHYTO == purpose, "TRVLCMIN"]
        ecdf_of_times = ecdf(times)
        correction = ecdf_of_times(1)
        weights_list$i = list(weights=ecdf_of_times, correction=correction)
        rm(times, ecdf_of_times, correction)
      }
    }
    print(paste("Weighting function obtained for", i))
  }
  rm(rest, n_city, n_state, purpose, i)
print("Constructing weights for walk-time distances to transit...")
print("(based on the AOI, where granularity is a function of NHTS samples occuring in the AOI)")
  tran = nhts_list_of_tables$just_transit
  n_city = sum(tran$HH_CBSA == city)
  n_state = sum(tran$HHSTATE == state)
  if(n_city >= ss){
    print("Weights for transit will be based on trips in the provided metro area")
    times = tran[tran$HH_CBSA == city, "TRACCTM"]
    ecdf_of_times = ecdf(times)
    correction = ecdf_of_times(1)
    weights_list$transit = list(ecdf_of_times, correction)
    rm(times, ecdf_of_times, correction)
  } else{
    if(n_state >= ss){
      print("Weights for transit will be based on trips in the provided metro area's state")
      times = tran[tran$HHSTATE == state, "TRACCTM"]
      ecdf_of_times = ecdf(times)
      correction = ecdf_of_times(1)
      weights_list$transit = list(ecdf_of_times, correction)
      rm(times, ecdf_of_times, correction)
    } else{
      print("Weights for transit will be based on trips in all available locations")
      times = tran[tran$HHSTATE == state, "TRACCTM"]
      ecdf_of_times = ecdf(times)
      correction = ecdf_of_times(1)
      weights_list$transit = list(weights=ecdf_of_times, correction=correction)
      rm(times, ecdf_of_times, correction)
    }
  }
  print("Weighting function obtained for transit")
  rm(tran, n_city, n_state, code, state)
print("Done! Weighting functions successsfully obtained for all non-transit indicators")
  e = Sys.time()
  difftime(e, s, units="mins")
  rm(e, s)
  return(weights_list)
}

rm(avg,ni,y)

y = ecdf(NYC_walkscores)
y()
y = ecdf(DC_walkscores)
y(50)
