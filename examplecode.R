#Example code

## compute 95% equal-tail credible interval without quantile()
kPostCDF <- cumsum(distribution)*gridwidth
lowerBound <- xgrid[which.min(abs(kPostCDF-0.025))]
upperBound <- xgrid[which.min(abs(kPostCDF-0.975))]
print(c(lowerBound,upperBound))

## normalzie posterior distribution, use exp if log unnormalized
normalized <- 1/gridwidth*exp(distribtion)/sum(exp(distribution))

## optim example 
toret <- optim(1,tempfunc, b = beta, errorval = 0, method = 'Brent', lower = 0, upper = 1,control=list(fnscale=-1))

# 95% equal-tail credible interval with quantile()
patientinterval <- quantile(patientdraws, probs=c(0.05,0.95))
