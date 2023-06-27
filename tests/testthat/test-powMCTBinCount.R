context("power calculation binary and count data")

## general options
mvt_control <- DoseFinding:::mvtnorm.control(maxpts=1e5, abseps = 0.0001)

## example (binary data)
candModList <- list(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1))

powA <- DoseFinding:::powMCTBinCount(rep(20,5), doses = c(0, 0.5, 1.5, 2.5, 4),
                       candModList=candModList, respModList=candModList,
                       placEffu = 0.1, maxEffu = 0.25, 
                       type = "binary_logit", option = "A",
                       alpha = 0.1, theta, control = mvt_control,
                       addArgs = list(scal = 4.8))
## externally calculated result in order Emax1,  Emax2,  sigEmax1, sigEmax2, betaMod
extern_pow <- c(0.6792, 0.7365, 0.825, 0.8007, 0.7139)
expect_equal(unname(powA), extern_pow, tolerance = 0.001)

## just test whether code fails (for testing option = "B" see further below)
powB <- DoseFinding:::powMCTBinCount(rep(20,5),doses = c(0, 0.5, 1.5, 2.5, 4),
                       candModList=candModList, respModList=candModList,
                       placEffu = 0.1, maxEffu = 0.25, 
                       type = "binary_logit", option = "B",
                       alpha = 0.1, theta, control = mvt_control,
                       addArgs = list(scal = 4.8))

## example (negative binomial data)
candModList <- list(emax = c(0.1, 0.5, 1, 2), sigEmax = rbind(c(1, 3), c(3, 3)),
                    betaMod = c(0.37, 0.74))
powA <- DoseFinding:::powMCTBinCount(c(100,50,50,100,100), doses = c(0, 0.5, 2, 4, 8), 
                       candModList=candModList, respModList=candModList,
                       placEffu = 0.6, maxEffu = -0.3, 
                       type = "negative_binomial", option = "A",
                       alpha = 0.025, theta = 1.25, control = mvt_control,
                       addArgs = list(scal = 9.6))
## externally calculated result in order Emax1,  Emax2,  Emax3,  Emax4, sigEmax1, sigEmax2, betaMod
extern_pow <- c(0.9035, 0.8816, 0.8691, 0.8518, 0.9264,  0.8518, 0.7913)
expect_equal(unname(powA), extern_pow, tolerance = 0.001)

## just test whether code fails (for testing option = "B" see further below)
powB <- DoseFinding:::powMCTBinCount(c(100,50,50,100,100), doses = c(0, 0.5, 2, 4, 8), 
                       candModList=candModList, respModList=candModList,
                       placEffu = 0.6, maxEffu = -0.3, 
                       type = "negative_binomial", option = "B",
                       alpha = 0.025, theta = 1.25, control = mvt_control,
                       addArgs = list(scal = 9.6))


## tests for option = "B"
## cannot validate against external results, so validate against manually calculated result
cVal <- qnorm(1-0.05)
n <- 50

## binary case
logit <- function(x)
  log(x/(1-x))
placEffu <- 0.4
maxEffu <- 0.3
p1 <- placEffu+maxEffu
p0 <- placEffu
nc_num <- logit(p1)-logit(p0)
nc_den <- sqrt(1/(n*p0*(1-p0))+1/(n*p1*(1-p1)))
delta <- nc_num/nc_den
power_bin <- 1-pnorm(cVal, delta, 1)
## compare against powMCTBinCount
candModList <- list(linear = NULL)
powB <- DoseFinding:::powMCTBinCount(c(n,n), doses = c(0, 1), 
                       candModList=candModList, respModList=candModList,
                       placEffu = placEffu, maxEffu = maxEffu, 
                       type = "binary_logit", option = "B",
                       alpha = 0.05, control = mvt_control,
                       addArgs = list(scal = 1))
expect_equal(unname(powB), power_bin, tolerance = 0.000001)

## negative binomial case
placEffu <- 1
maxEffu <- -0.5
theta <- 2
r0 <- placEffu
r1 <- placEffu + maxEffu
nc_num <- -1*(log(r1)-log(r0)) # multiply with -1 (decreasing)
v <- (theta+c(r0, r1))/(theta*c(r0, r1))
nc_den <- sqrt(sum(v/c(n,n)))
delta <- nc_num/nc_den
power_nb <- 1-pnorm(cVal, delta, 1)

powB <- DoseFinding:::powMCTBinCount(c(n,n), doses = c(0, 1), 
                       candModList=candModList, respModList=candModList,
                       placEffu = placEffu, maxEffu = maxEffu, 
                       type = "negative_binomial", option = "B",
                       alpha = 0.05, theta = theta, control = mvt_control,
                       addArgs = list(scal = 1))
expect_equal(unname(powB), power_nb, tolerance = 0.000001)

## tests for contrast matrix handed over (just test whether code fails)
contMat <- rbind(rep(-1, 4), diag(4))
respModList <- list(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1))
powA <- DoseFinding:::powMCTBinCount(n=rep(20,5), doses = c(0, 0.5, 1.5, 2.5, 4),
                       candModList=NULL,
                       respModList=respModList,
                       placEffu = 0.1, maxEffu = 0.25, 
                       type = "binary_logit", 
                       alpha = 0.1, theta, control = mvt_control,
                       contMat = contMat,
                       addArgs = list(scal = 4.8))

respModList <- list(emax = c(0.1, 0.5, 1, 2), sigEmax = rbind(c(1, 3), c(3, 3)),
                    betaMod = c(0.37, 0.74))
contMat <- rbind(rep(1, 4), -diag(4))
powA <- DoseFinding:::powMCTBinCount(c(100,50,50,100,100), doses = c(0, 0.5, 2, 4, 8), 
                       respModList=respModList,
                       placEffu = 0.6, maxEffu = -0.3, 
                       type = "negative_binomial", 
                       alpha = 0.025, theta = 1.25, control = mvt_control,
                       contMat = contMat,
                       addArgs = list(scal = 9.6))
