set.seed(12345)
###### Q1
##b)

x <- 33
gridstep <- 0.001
thetaGrid <- seq(0,1,gridstep)
unnormpost <- dbinom(x,50,thetagrid)*(thetaGrid>0.3)*(thetaGrid>0.7)
postB <- 1/gridstep*unnormpost/sum(unnormpost)
plot(thetaGrid,postB, type = 'l')
adens <- dbeta(thetaGrid,x+1,51-x)
lines(thetaGrid, adens,type = 'l', col = 2)
probA <- pbeta(0.5,x+1,51-x)
probB <- sum(postB[thetaGrid<=0.5]*gridstep)
print(c(probA,probB))
