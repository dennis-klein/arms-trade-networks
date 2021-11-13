library(RSiena)


### LOAD DATA

load("data/out/sipri_tiv.RData") # sipri_tiv: networks on arms trade amount T=69
load("data/out/atop_alliance.RData") # atop_alliance: indicates whether countries are allies
load("data/out/pop.RData") # Population
load("data/out/real_gdp_p_c.RData") # GDP p.c.
load("data/out/cepii_trade.RData")
# mil exp

### PREPARE DATA
N.countries <- 257
T.periods <- 69
#N.countries <- dim(sipri_tiv[[1]])[1]
#T.periods <- length(sipri_tiv)

# convert TIV amount to trade dummy
# later we will try to model flow data
armstrade.net <- array(0, dim = c(N.countries, N.countries, T.periods))
for (t in 1:T.periods) {
  armstrade.net[,,t] <- ifelse(sipri_tiv[[t]] > 0, 1, 0)
}
mode(armstrade.net) <- "integer"

# replace zero values with NAN in Pop, GDP p.c.
pop[pop == 0] <- NA
real_gdp_p_c[real_gdp_p_c == 0] <- NA

# log transformation for Pop and GDP p.c.
# also helps to scale down data in the recommended standard deviation area of 0.1 to 10
# TODO: check whether log scale OK here
log_pop <- log(pop+1)
log_gdppc <- log(real_gdp_p_c+1)

# cepii trade data
trade <- array(0, dim = c(N.countries, N.countries, T.periods))
for (t in 1:T.periods) {
  trade[,,t] <- cepii_trade[[t]]
}
# cut last period (model requirement of T.periods - 1 observations)
trade <- trade[,,1:T.periods-1]
# log transformation
log_trade <- log(trade+1)

# NOTE: maybe consider change in trade


### BUILD THE SAOM MODEL
# define dependent variables
model.dep <- sienaDependent(armstrade.net)

# define actor covariates (coCovars and varCovars)
model.varc.log_pop <- varCovar(log_pop)
model.varc.log_gdppc <- varCovar(log_gdppc)

# define dyadic covariates (coDyadCovar and varDyadCovar)
model.vardyadc.log_trade <- varDyadCovar(log_trade)

# SIENA data object
model.data <- sienaDataCreate(
  model.dep,
  model.varc.log_gdppc
  #model.varc.log_pop,
  #model.vardyadc.log_trade
)

#print01Report(model.data)

# create SIENA effects
model.effects <- getEffects(model.data)

# create SIENA algorithm
model.algo <- sienaAlgorithmCreate(projname = "SAOM_first_try")

# estimate SAOM model
model.est <- siena07(model.algo, data = model.data, effects = model.effects)




