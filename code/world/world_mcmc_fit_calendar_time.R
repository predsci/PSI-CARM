rm(list = ls())

library(squire)
library(countrycode)
library(EpiEstim)
library(incidence)
library(ggplot2)
library(grid)
library(gridExtra)
library(coda)

## 
## Provide start/end index for location fitting if only start is provided end will be set to last one in dataset
##

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    nc_start = 1
    nc_end = NULL
} else if (length(args) == 1) {
    nc_start = as.numeric(args[1])
    nc_end = NULL
} else if (length(args) == 2) {
    nc_start = as.numeric(args[1])
    nc_end   = as.numeric(args[2])
} else {
    stop("Only 0, 1 or 2 arguments must be supplied\n", call=FALSE)

}

topn = 120
## 
## social mixing: either uniform or not 
##

uniform = FALSE

##
## Set the last day for modeling 
##

date_last = as.Date("2020-04-15") + 75
print("Fitting data until:")
print(date_last)

days_from_first_death = date_last
##
## Number of unique values for R(t) - loop over all these values- will need to increase this with date
##

## for 04-15 use 2:3, for the next two use 2:4 and for the rest 2:6

if (date_last == as.Date("2020-04-15")) {
	nb_vec = 2:3
} else if (date_last >  as.Date("2020-04-15") & date_last <= (as.Date("2020-04-15")+30)) {
	nb_vec = 2:4
} else {
	nb_vec = 2:5
}

##
##	Length of MCMC chain
##
	
mcmclen = 1e6

##
## Create sub-directories for results
## 
dir_names_data_files = paste0('/fit_files_',nb_vec,'/')
dir_names_plot_files = paste0('/plots_',nb_vec,'/')

dir_name_days_from_first_death = paste0('fits_from_death1_to_',days_from_first_death,'/')

workdir = getwd()

if(!dir.exists(dir_name_days_from_first_death)) dir.create(dir_name_days_from_first_death,'/')

outdir_vec  = paste0(workdir,'/',dir_name_days_from_first_death, dir_names_data_files)
plotdir_vec = paste0(workdir,'/',dir_name_days_from_first_death, dir_names_plot_files)


for (jj in 1:length(nb_vec)) {
	if(!dir.exists(outdir_vec[jj])) dir.create(outdir_vec[jj])
	if(!dir.exists(plotdir_vec[jj])) dir.create(plotdir_vec[jj])
}


##
## Load the Fortran codes
##

dyn.load("../src/tanhrt.so")

dyn.load("../src/mcmc.so")

dyn.load("../src/detcovid.so")

dyn.load("../src/stochcovid.so")

## 
## Test that Loading worked
## 

if(!is.loaded("tanhrt")) {
  cat("\nFailed to load ../src/tanhrt.so\n Code will exit \n\n" )
  quit()
}


if(!is.loaded("mcmc")) {
  cat("\nFailed to load ../src/mcmc.so\n Code will exit \n\n" )
  quit()
}


if(!is.loaded("detcovid")) {
  cat("\nFailed to load ../src/detcovid.so\n Code will exit \n\n" )
  quit()
}


if(!is.loaded("stochcovid")) {
  cat("\nFailed to load ../src/stochcovid.so\n Code will exit \n\n" )
  quit()
}


trim.data.out <- function(longvec, tiny=1e-10) {
	noobs <- length(longvec)
	first <- 1
	while (longvec[first] < tiny) {
		first = first + 1
	}
	first
}

##
## Logical flags for plotting
##

plot = TRUE
plot_wt = FALSE
save_analysis_csv = FALSE 

##
## Load data 
## 
DRAFT_db_file = "../../data/World_data.rds"

out_list = readRDS(file=DRAFT_db_file)

# Remove United States states Mariana Islands from database structure

remove_index = which(out_list$locs$Country.Region == "United States" & out_list$locs$Province.State == "Northern Mariana Islands")
out_list$locs = out_list$locs[-remove_index, ]
out_list$time_series = out_list$time_series[-remove_index]
out_list$pop_info = out_list$pop_info[-remove_index]

dates <- out_list$time_series[[1]]$date

## retrieve cases and deaths from out_list

cases <- as.data.frame(matrix(data=0, nrow=nrow(out_list$locs), ncol=length(dates)))
names(cases) = format(dates, "%Y-%m-%d")
for (ii in 1:nrow(cases)) cases[ii, ] = out_list$time_series[[ii]]$cumul_cases

deaths <- as.data.frame(matrix(data=0, nrow=nrow(out_list$locs), ncol=length(dates)))
names(deaths) = format(dates, "%Y-%m-%d")
for (ii in 1:nrow(deaths)) deaths[ii, ] = out_list$time_series[[ii]]$cumul_deaths


remove_index <- which(dates > date_last)

if (length(remove_index) != 0) {
        dates <- dates[-remove_index]
        cases <- cases[, -remove_index]
        deaths <- deaths[, -remove_index]
}

##
## Get the country, province/state, population and age distributions
##

country <- out_list$locs$Country.Region
prov_state <- out_list$locs$Province.State
population <- out_list$locs$population
age_data <- out_list$pop_info

ndates = length(dates)


## Now reorder and take only topn location

iorder <- order(deaths[,ndates], decreasing = TRUE)

country <- country[iorder]
prov_state <- prov_state[iorder]
population <- population[iorder]
age_data <- age_data[iorder]

cases <- cases[iorder, ]
deaths <- deaths[iorder, ]

# now take only topn

cases <- cases[1:topn, ]
deaths <- deaths[1:topn, ]

country <- country[1:topn]
prov_state <- prov_state[1:topn]
population <- population[1:topn]
age_data <- age_data[1:topn]

## we choose to draw:

config <- make_config(list(mean_si = 6.48, std_mean_si = 3.83, min_mean_si = 2.48, max_mean_si = 10.48, std_si = 10, std_std_si = 1, min_std_si = 1, max_std_si = 19))

nc = dim(cases)[1]

cat("Total Number of locations: ", nc, "\n")

nReals = 1000

I0total = 10

dt=1
	
vecPS = c(1.61e-05, 6.95e-05, 0.000309, 0.000844, 0.00161, 0.00595, 0.0193, 0.0428, 0.078)
vecPM = c(0, 0.00041, 0.014, 0.0343, 0.0425, 0.0816, 0.118, 0.166, 0.184)

if (is.null(nc_end)) nc_end = nc

##
## Sanity check to ensure we are running on the correct range of states
##
if (nc_end > nc) nc_end = nc

if (nc_start > nc_end ) {
    stop(paste0("nc_start = ", nc_start, " is gt than nc_end = ", nc_end, "\n"), call=FALSE)
}

##
## Prepare the analysis table - this will also need to be updated after all the countries are done 
##

cat("\n Processing Locations from ", nc_start, " to ", nc_end, "\n")

for (ii in nc_start:nc_end) {

	mycountry = country[ii]
	
	mystate = prov_state[ii]

	myname <- paste0(mycountry, " ", mystate)

	# retrieve cum. death for this country
	
	obs_all = deaths[ii, ]

	# incidence death
	inc_all = diff(as.numeric(deaths[ii, ]))
	
	inc_all = c(0, inc_all)
	
	cat(ii, mycountry, mystate, "\n")
	
	myiso3 = countrycode(sourcevar = mycountry, destination = "iso3c", origin = "country.name")
	
	age_dist = age_data[[ii]]$age_dist

	nAges = length(age_dist)

	vecPA <- rep(0.4, nAges)

	cat("\n Processing Location Number", ii, "Name: ", myname, "Population: ", population[ii],  "\n\n")


	# find the index and date of first death
	ind.dd1 <- which(obs_all > 0)[1]
	dd1 <- dates[ind.dd1]

	# use Corry procedure to get R(t)
	daily <- diff(as.numeric(cases[ii, ]))
	daily <- c(0, daily)
	first <- trim.data.out(longvec = daily, tiny = 10)

	dates.trim <- dates[first:ndates]
	daily.trim <- daily[first:ndates]

	for (jj in 1:length(daily.trim)) daily.trim[jj] = max(daily.trim[jj], 0)

	mydata <- as.incidence(daily.trim, dates = dates.trim)
	res_parametric_si <- estimate_R(mydata, method = "uncertain_si", config = config)

	wt.dates <- res_parametric_si$dates
	rt <- res_parametric_si$R$Mean
	t_start <- res_parametric_si$R$t_start
	t_end <- res_parametric_si$R$t_end

	week = 7
	x = (week + 1):length(wt.dates)
	y = rt
	lo = loess(y ~ x, control = loess.control(surface = "direct"), model= TRUE)
	
	x.expand = 1:length(wt.dates)
	# the length of 'srt' is the same as dates.trim
	srt = predict(lo, newdata = x.expand)
	# Do not let it go down
	if (srt[1] < srt[week + 1]) {
		srt[1:week] = srt[week + 1]
	}

	# need to fill back dates with initial value
	srt.all = rep(srt[1], length(dates))
	srt.all[first:length(dates)] <- srt

	if (plot_wt == TRUE) {
		filename = paste0(workdir,"/plots/wt_", myname, ".pdf")
		pdf(file = filename)
		plot(res_parametric_si, legend = FALSE, options_I = list(col = "coral"), options_R = list(col = "blue"), options_SI = list(col = "grey", transp = 0.2))
		dev.off()
	}


	## now that we have an estimat for R0 we can see how long does it take to get to a cumulative number of 1,000 cases.
	list_rt <- srt.all
	list_dates <- dates
	list_dates_trim <- dates.trim
	list_rt_trim <- srt
	
	cat(myname, srt, "\n")

	Ntotal = population[ii]

	nAges = length(age_dist)
	vecPA <- rep(0.4, nAges)
	vecN = round(Ntotal * age_dist)
	vecI0 <- round(I0total * age_dist)
	srt <- list_rt
	R0 = srt[1]


	ndays = 90
	xvals = seq(0, ndays, 1)
	vecTcalc = seq(0, ndays, dt)
	vecRtrel = rep(1, length(vecTcalc)) #srt.all[(ndates-ndays):ndates]

	trig_pres = 99999
	trig_day = 99999
	icu_cap = 0.9999
	trickle = 0
	D_E=2
	D_I1=3
	D_I2=3
	D_HR=20
	D_HD=23
	D_ICU=rep(7, nAges)
	sevBchange=0 # 0 = FALSE 1 = TRUE
	vecInfNess=rep(1, nAges)
	vecSusTy=rep(1, nAges)

	
	if (uniform) {
		matCt = matrix(1/nAges, nrow = nAges, ncol = nAges)

		matCtClosure = matrix(1/nAges, nrow = nAges, ncol = nAges)
	} else {

		mata = get_mixing_matrix(iso3c = myiso3)
		sqr_pop = get_population(iso3c = myiso3)$n
		pop_copy = sqr_pop
		n = dim(mata)[1]
		nsquire = length(sqr_pop)
		if (n < nsquire) 
			sqr_pop[n] = sum(sqr_pop[n:nsquire])
		sqr_pop = sqr_pop[-nsquire]
		nsquire = length(sqr_pop)
		mats = mata * 0

		for (jj in 1:nsquire) {

			for (kk in 1:nsquire) {
				mats[jj, kk] = 1/(2 * sqr_pop[jj]) * (mata[jj, kk] * sqr_pop[jj] + mata[kk, jj] * sqr_pop[kk])
			}

		}

		mat_psi = array(0, c(nAges, nAges))
		icount = 1
		mat_tmp = array(0, c(nAges, nsquire))
		for (jj in 1:(nAges - 2)) {
			mat_tmp[jj, 1:nsquire] = (mats[icount, 1:nsquire] * pop_copy[icount] + mats[(icount + 1), 1:nsquire] * pop_copy[(icount + 1)])/(pop_copy[icount] + pop_copy[(icount + 
				1)])
			icount = icount + 2
		}

		for (jj in (nAges - 1):nAges) {
			mat_tmp[jj, ] = mats[icount, ]
			icount = icount + 1
		}

		icount = 1
		for (jj in 1:(nAges - 2)) {
			mat_psi[1:nAges, jj] = (mat_tmp[1:nAges, icount] * pop_copy[icount] + mats[1:nAges, (icount + 1)] * pop_copy[(icount + 1)])/(pop_copy[icount] + pop_copy[(icount + 
				1)])
			icount = icount + 2
		}
		for (jj in (nAges - 1):nAges) {
			mat_psi[, jj] = mat_tmp[, icount]
			icount = icount + 1
		}

		#vecN - psi age distribution, need to symmetrize the matrix! 
		pop_psi = age_data[[ii]]$age_dist

		mix_mat = mat_psi * 0
		for (jj in 1:nAges) {

			for (kk in 1:nAges) {
				mix_mat[jj, kk] = 1/(2 * vecN[jj]) * (mat_psi[jj, kk] * vecN[jj] + mat_psi[kk, jj] * vecN[kk])
			}

		}
		matCt = mix_mat
		matCtClosure = mix_mat
		
		# In case of an NA - resort to uniform mixing and print a warning 
		if(any(is.na(mix_mat))) {
			cat('\n\n',"**** Warning: Found NAs in Mixing Matrix, Using Uniform Mixing Instead ****", '\n\n')
			matCt = matrix(1/nAges, nrow = nAges, ncol = nAges)
			matCtClosure = matrix(1/nAges, nrow = nAges, ncol = nAges)			
		}
	}

	scLim = c(99999,99999)

	nTimes = length(vecTcalc)

	rtn_inf <- array(0, dim=c(nTimes,nAges, nReals))
	rtn_icu <- array(0, dim=c(nTimes,nAges,nReals))
	rtn_ded_cum <- array(0.0, dim=c(nTimes,nAges,nReals))  

	daily_inf <- daily_icu <-  array(0, dim = c(ndays, nAges, nReals))
	daily_cded <-  array(0.0, dim = c(ndays, nAges, nReals))
	rt_daily <- array(0.0, ndays)
		
	Rval = R0
	tn = 0.0
	new.seed <- as.integer(runif(1)*2e4)
	
	# Find day 1 - see when we have more than 1,000 cases 
	
	y1 <- .Fortran('stochcovid', nb = as.integer(1), Rval = as.double(Rval), tn = as.double(tn), vecN = as.integer(vecN), vecI0 = as.integer(vecI0), vecInfNess = as.double(vecInfNess), vecSusTy = as.double(vecSusTy), vecPS = as.double(vecPS), vecPM = as.double(vecPM), vecPA = as.double(vecPA), matCt = as.double(matCt), matCtClosure=as.double(matCtClosure), scLim = as.double(scLim), vecTcalc = as.double(vecTcalc), vecRtrel = as.double(vecRtrel), D_E = as.double(D_E), D_I1 = as.double(D_I1), D_I2 = as.double(D_I2), D_HR = as.double(D_HR), D_HD = as.double(D_HD), D_ICU = as.double(D_ICU), trig_pres = as.integer(trig_pres), icu_cap = as.double(icu_cap), trickle = as.integer(trickle), sevBchange = as.integer(sevBchange), nAges = as.integer(nAges), nTimes = as.integer(nTimes), rtn_inf = as.integer(rtn_inf), rtn_icu = as.integer(rtn_icu), rtn_ded_cum = as.double(rtn_ded_cum), nReals = as.integer(nReals), iseed = as.integer(new.seed), ndays = as.integer(ndays), daily_inf = as.integer(daily_inf),  daily_icu = as.integer(daily_icu),  daily_cded = as.double(daily_cded), rt_daily = as.double(rt_daily))

   
	rtn_inf = array(y1$daily_inf, c(ndays, nAges, nReals))
	
	
	cdays <- rep(NA, nReals)
	for (jj in 1:nReals) {
		cdays[jj] <- which(rowSums(rtn_inf[, ,jj]) > 1000)[1]

	}

	if (all(is.na(cdays))) {
		first_case_day <- which(daily > 0)[1]
		day1 = max(1, first_case_day)
		save_day1 = day1
		day1.date <- dates[day1]
	} else {
		day1 <- ind.dd1 - round(median(cdays, na.rm = TRUE))
		day1 = max(1, day1)
		save_day1 = day1
		day1.date <- dates[day1]
	}

	myname = mycountry
	if (mystate != "" & mystate != mycountry) 
		myname = paste(mycountry, ", ", mystate)
		
	country_name = myname
	state_name = mystate
	country_pop = Ntotal
	country_age = age_dist
	country_day1 = as.Date(day1.date, format = "%Y-%m-%d")

	country_all_death_dates = dates
	country_all_death = obs_all
	country_all_cases = cases[ii, ]
	country_all_inc_death = inc_all
	
	country_iso3 = myiso3
	
	# Now the fit - and for something new - we are going to loop over all the values of nb_vec
	
	for (inb in 1:length(nb_vec)) {
		
		nb= nb_vec[inb]
		outdir = outdir_vec[inb]
		plotdir = plotdir_vec[inb]
		
		cat("\n Starting MCMC Procedure for nb = ", nb, '\n\n')
	
	I0total = 10
	vecN = round(Ntotal * age_dist)
	vecI0 <- round(I0total * age_dist)

	## Here change to model only a subset of the data 
	day_first_death <- which(obs_all >= 1)[1]
	
	obs <- as.numeric(obs_all[save_day1:length(as.numeric(obs_all))])

	obs_cases <- as.numeric(cases[ii, save_day1:length(as.numeric(obs_all))])

	obs_inc <- as.numeric(inc_all[save_day1:length(as.numeric(obs_all))])
		
	gama_obs_inc = lgamma((obs_inc + 1))
	
	# take care of infinity values which the weights will zero later
	
	gama_obs_inc[is.infinite(gama_obs_inc)] <- 0
	
	## weights 
	
	wght = rep(1, length(obs_inc))
	
	ind_0 <- which(obs_inc <= 0)
	
	wght[ind_0] = 0.0
	
	##
	
	day1 = country_day1

	dates2_model <- subset(country_all_death_dates, country_all_death_dates >= day1)

	dates2_model <- dates2_model[1:length(obs_inc)]
	
	ndays <- length(dates2_model)

		Rval = rep(0, nb)
		tn = rep(0, nb)

		knots = nb + 1
		step = floor(length(rt)/knots)
		for (i in 1:nb) {
			start = 1 + (i - 1) * step
			end = step * i

			Rval[i] = mean(rt[start:end])
		}

		step = floor(length(obs)/nb)
		tn[1:nb] = round(runif(nb, step - 2, step + 2))
    
    Rval_min = rep(0.5, nb)
	Rval_max = rep(5.0, nb)
	Rval_max[2:nb] = 2.0
	tn_min = rep(7, nb)
	
	tn_min[nb] = 7
	
	tn_max = rep(ndays-7, nb)

	# Sanity check - make sure initial guess is within min/max values
	
	for (i in 1:nb) {
		if(Rval[i] >= Rval_max[i] || Rval[i] <= Rval_min[i]) Rval[i] = runif(1, Rval_min[i], Rval_max[i])
		if (tn[i] >= tn_max[i] || tn[i] <= tn_min[i]) tn[i] = runif(1, tn_min[i], tn_max[i])
	}
	
	tn[nb] = 0
    drval = rep(0.4, nb)
    drval[1] = 0.5
    dtn = rep(1.0, nb)
    
	nparam = nb * 2

	cat("\n MCMC Chain Length: ", mcmclen, "\n\n")
	cat("\n Fitting ",ndays," days \n\n")
	
	# Number of trials to keep
	
    nlines = mcmclen # 1e4
	theta_mcmc = array(0, c(nlines, nparam))
	
	ndays = length(obs)

	xvals = seq(0, ndays, 1)
	vecTcalc = seq(0, ndays, dt)
	vecRtrel = rep(1, length(vecTcalc)) 
	nTimes = length(vecTcalc)
    
    new.seed <- as.integer(runif(1)*2e4)

	out <- .Fortran('mcmc', nb = as.integer(nb), Rval= as.double(Rval), tn = as.double(tn), Rval_max = as.double(Rval_max), Rval_min= as.double(Rval_min), tn_max=as.double(tn_max), tn_min = as.double(tn_min), Rval_step = as.double(drval), tn_step = as.double(dtn), vecN = as.double(vecN), vecI0 = as.double(vecI0), vecInfNess = as.double(vecInfNess), vecSusTy = as.double(vecSusTy), vecPS = as.double(vecPS), vecPM = as.double(vecPM), vecPA = as.double(vecPA), matCt = as.double(matCt), matCtClosure=as.double(matCtClosure), scLim = as.double(scLim), vecTcalc = as.double(vecTcalc), vecRtrel = as.double(vecRtrel), D_E = as.double(D_E), D_I1 = as.double(D_I1), D_I2 = as.double(D_I2), D_HR = as.double(D_HR), D_HD = as.double(D_HD), D_ICU = as.double(D_ICU), trig_pres = as.integer(trig_pres), icu_cap = as.double(icu_cap), trickle = as.double(trickle), sevBchange = as.integer(sevBchange), nAges = as.integer(nAges), nTimes = as.integer(nTimes), obs_inc = as.double(obs_inc), gama_obs_inc = as.double(gama_obs_inc), wght = as.double(wght), ndays = as.integer(ndays),theta = as.single(theta_mcmc), mcmclen = as.integer(mcmclen), nlines = as.integer(nlines), iseed = as.integer(new.seed))
	

	theta_mcmc = array(out$theta, c(nlines, (2 * nb)))

	# convert RSS to AIcc - note since our RSS already divided by the number of data points we do not divide again here
	# Formula is AIC = 2 * k + n * ln(RSS/n)
	# RSS = sum_over_i(data_i-model_i)^2
	#

	## Since we are using the possion now (-ln(L))
	
	## AICc = 2*nopt - 2*ln(L) + (2*nopt*(nopt+1))/(ndays-nopt-1)
	nopt = nb + nb - 1

	theta_mcmc[, (2 * nb)] = 2 * nopt + 2 * theta_mcmc[, (2 * nb)] + (2 * nopt * (nopt + 1))/(ndays - nopt - 1)


	Rval_mcmc = theta_mcmc[, 1:nb]
	tn_mcmc = theta_mcmc[, (nb + 1):(2 * nb - 1)]

	iburn = 1/5 * mcmclen
	ithin = mcmclen/nlines
	
	colnames(theta_mcmc) <- c(paste0("Rval", 1:nb), paste0('tn_', 1:(nb-1)), 'AICc')

	results = mcmc(data = theta_mcmc, start = 1, end = mcmclen, thin = ithin)
	
	country_mcmc_summary = summary(results[,(nb*2)], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE))
	
	print(summary(results[,(nb*2)], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
	
	ess <- effectiveSize(results[,1:(nparam-1)])
	
	print("ESS")
	print(ess)
	##
	## These are the best parameters
	##

	ll = theta_mcmc[,(2*nb)]
	
	ibest = which.min(ll)

	cat("\n Theta Best: ", theta_mcmc[ibest, ], "\n")
			
	Rval = rep(0, nb)
	tn = rep(0, nb)
	
	Rval[1:nb] = theta_mcmc[ibest, 1:nb]
	tn[1:(nb-1)] = theta_mcmc[ibest, (nb+1):(2*nb-1)]
	tn[nb] = 0.0
		
    rtn_inf = rtn_icu = rtn_ded_cum = array(0.0, c(nTimes, nAges))
	
    daily_inf <- daily_icu <- daily_cded <-  array(0.0, dim = c(ndays, nAges))
	
    rt_daily <- array(0, ndays)

	y2 <- .Fortran('detcovid', nb = as.integer(nb), Rval = as.double(Rval), tn = as.double(tn), vecN = as.double(vecN), vecI0 = as.double(vecI0), vecInfNess = as.double(vecInfNess), vecSusTy = as.double(vecSusTy), vecPS = as.double(vecPS), vecPM = as.double(vecPM), vecPA = as.double(vecPA), matCt = as.double(matCt), matCtClosure=as.double(matCtClosure), scLim = as.double(scLim), vecTcalc = as.double(vecTcalc), vecRtrel = as.double(vecRtrel), D_E = as.double(D_E), D_I1 = as.double(D_I1), D_I2 = as.double(D_I2), D_HR = as.double(D_HR), D_HD = as.double(D_HD), D_ICU = as.double(D_ICU), trig_pres = as.double(trig_pres), icu_cap = as.double(icu_cap), trickle = as.double(trickle), sevBchange = as.integer(sevBchange), nAges = as.integer(nAges), nTimes = as.integer(nTimes), rtn_inf = as.double(rtn_inf), rtn_icu = as.double(rtn_icu), rtn_ded_cum = as.double(rtn_ded_cum), ndays = as.integer(ndays), daily_inf = as.double(daily_inf), daily_icu = as.double(daily_icu), daily_cded = as.double(daily_cded), rt_daily = as.double(rt_daily))


	ddaily_inf  = array(y2$daily_inf , c(ndays, nAges))
	ddaily_icu  = array(y2$daily_icu , c(ndays, nAges))
	ddaily_cded = array(y2$daily_cded, c(ndays, nAges))

	rt_daily = y2$rt_daily
		
	#cat('\n rt_daily', rt_daily, '\n\n')

	# Update this after the fit including day2_start
	country_llk = theta_mcmc[ibest,(2*nb)]
	country_rt = round(rt_daily, digits = 2)
	country_Rval = Rval
	country_tn = tn

	country_ess = ess
	
	nReals = 100

	rtn_inf <- array(0, dim=c(nTimes,nAges, nReals))
	rtn_icu <- array(0, dim=c(nTimes,nAges,nReals))
	rtn_ded_cum <- array(0.0, dim=c(nTimes,nAges,nReals))  

	daily_inf <- daily_icu <-  array(0, dim = c(ndays, nAges, nReals))
	daily_cded <-  array(0.0, dim = c(ndays, nAges, nReals))
		
	new.seed <- as.integer(runif(1)*2e4)
	y1 <- .Fortran('stochcovid', nb = as.integer(nb), Rval = as.double(Rval),tn = as.double(tn), vecN = as.integer(vecN), vecI0 = as.integer(vecI0), vecInfNess = as.double(vecInfNess), vecSusTy = as.double(vecSusTy), vecPS = as.double(vecPS), vecPM = as.double(vecPM), vecPA = as.double(vecPA), matCt = as.double(matCt), matCtClosure=as.double(matCtClosure), scLim = as.double(scLim), vecTcalc = as.double(vecTcalc), vecRtrel = as.double(vecRtrel), D_E = as.double(D_E), D_I1 = as.double(D_I1), D_I2 = as.double(D_I2), D_HR = as.double(D_HR), D_HD = as.double(D_HD), D_ICU = as.double(D_ICU), trig_pres = as.integer(trig_pres), icu_cap = as.double(icu_cap), trickle = as.integer(trickle), sevBchange = as.integer(sevBchange), nAges = as.integer(nAges), nTimes = as.integer(nTimes), rtn_inf = as.integer(rtn_inf), rtn_icu = as.integer(rtn_icu), rtn_ded_cum = as.double(rtn_ded_cum), nReals = as.integer(nReals), iseed = as.integer(new.seed), ndays = as.integer(ndays), daily_inf = as.integer(daily_inf),  daily_icu = as.integer(daily_icu),  daily_cded = as.double(daily_cded),  rt_daily = as.double(rt_daily))
	
	sdaily_inf  = array(y1$daily_inf , c(ndays, nAges, nReals))
	sdaily_icu  = array(y1$daily_icu , c(ndays, nAges, nReals))
	sdaily_cded = array(y1$daily_cded, c(ndays, nAges, nReals))

	rtn_ded_cum=array(y1$rtn_ded_cum, c(nTimes, nAges, nReals))
	
	## Calculate the median 
	
	med_cded = med_inf = array(0, c(ndays, nAges))
	
	for (i in 1:ndays) {
		for (j in 1:nAges) {
			med_cded[i,j] = mean(sdaily_cded[i,j,])
            med_inf[i,j]  = mean(sdaily_inf[i,j,])
		}
	}	

	## sum over ages 
	
	model_cded  = array(NA, c(ndays, nReals))
	for (i in 1:ndays) {
		for (j in 1:nReals) {
			model_cded[i, j] = sum(sdaily_cded[i,, j])
		}
	}

	string_ii = formatC(ii, flag=0, width=3)
	
	if (plot == TRUE) {
		file1 = paste0(plotdir,"/fits_", string_ii, ".pdf")
		pdf(file = file1)
		if (nb == 1) {
			par(mfrow = c(2, 1))
		} else if (nb == 2 || nb == 3) {
			par(mfrow = c(2, 2))
		} else if (nb == 3 || 4) {
			par(mfrow = c(3, 2))
		} else {
			par(mfrow = c(4, 2))
		}

		cat("\n\n Creating Plot: ", file1, "\n\n")
		
		pl_obs_inc = obs_inc
		pl_obs_inc[obs_inc < 0] <- NA
		plot(dates2_model, pl_obs_inc, type = "p", col = "red", ylab = "Deaths", main = myname, bty = "n", xlab = 'Time')
		for (i in seq(from=1, to = nReals, length = 100)) lines(dates2_model, c(0,diff(rowSums(sdaily_cded[, , i]))), col = "grey")
		lines(dates2_model, c(0, diff(rowSums(med_cded))), col = "blue")
		lines(dates2_model, pl_obs_inc, type = "p", col = "red") 
		par(new = TRUE)
		plot(dates2_model, rt_daily, col = "black", lwd = 2, xaxt = "n", yaxt = "n", xlab ='', ylab = '', type = 'l', ylim = c(0.5, 5))
		axis(4)	
			
		plot(dates2_model, obs, type = "p", col = "red", ylab = "Cumulative Deaths", main = myname, bty = "n", xlab = 'Time')
		for (i in seq(from=1, to = nReals, length = 100)) lines(dates2_model, rowSums(sdaily_cded[, , i]), col = "grey")
		lines(dates2_model, obs, type = "p", col = "red",cex = 1.5)
		lines(dates2_model, rowSums(med_cded), col = "blue")	
		par(new = TRUE)
		plot(dates2_model, rt_daily, col = "black", lwd = 2, xaxt = "n", yaxt = "n", xlab ='', ylab = '', type = 'l', ylim = c(0.5, 5))
		axis(4)				
	
		for (i in 1:nb) hist(Rval_mcmc[, i], col = "coral", main = paste0("Rval Number: ", i), xlab = 'R_i')

				
		dev.off()

	}
 
	# Save all the results for this location   

	country_epi_data = list(country_name = country_name, state_name = state_name, country_iso3 = country_iso3, country_pop = country_pop, country_age = country_age, country_day1 = country_day1, country_rt = country_rt, country_ess = country_ess, country_model_cded = model_cded, 
	country_Rval = country_Rval, country_tn = country_tn, country_all_death_dates = country_all_death_dates, country_all_death = country_all_death, 
	country_all_cases = country_all_cases, country_all_inc_death  = country_all_inc_death, country_mcmc_summary = country_mcmc_summary, country_llk = country_llk)

	file2 = paste0(outdir,"/country_epi_data_",string_ii,".rds")
	saveRDS(country_epi_data, file = file2)
	cat("Saved Epi Data to: ", file2, "\n")

	file4 = paste0(outdir,"/mcmc_",string_ii,".rds")
	saveRDS(results, file = file4)
	cat("Saved MCMC Data to: ", file4, "\n")
	} # end of loop over nb_vec
	

} ## End of loop over locations

##
##
