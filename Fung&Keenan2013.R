# pmdSamplingDistYiN function conversion to R code

# Fung & Keenan 2013

# 1

pmfSamplingDistYiN <- function(M, NN, p_i, P_ii, yiN){
#   M <- 1000
#   NN <- 30
#   p_i <- 0.2
#   P_ii <- p_i^2
#   yiN <- p_i * (2*NN)
  
  MaxFunc <- max((yiN/2) - (M*p_i) + (M*P_ii), yiN - NN, 0)
  MinFunc <- min(M*P_ii, yiN/2, M-NN+yiN-(2*M*p_i)+(M*P_ii))
  LowerBound <- ceiling(MaxFunc)
  UpperBound <- floor(MinFunc)
  
  # Calculate P(YiN = yiN) according to eqn 9
  Numerator1 <- 0
  
  for(i in LowerBound:UpperBound){
    Summand <- choose(M*P_ii, i)*choose(2*M*(p_i-P_ii), yiN-(2*i))*
      choose(M+(M*P_ii)-(2*M*p_i), NN+i-yiN)
    Numerator1 <- Numerator1 + Summand
  }
  
  probOut <- Numerator1/choose(M, NN)
  return(probOut)
}

# 

# 2

AcceptanceRegion <- function(M, NN, p_i0, P_ii0, yiNobs, alpha){
  
  # find the lower bound of the acceptance region
  SumProb <- 0
  yiNlow <- 0
  dummy1 <- 0
  while(SumProb <= (alpha/2)){
    dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNlow)
    SumProb <- SumProb + dummy1
    yiNlow <- yiNlow + 1
  }
  yiNlow <- yiNlow - 1
  dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNlow)
  SumProb <- SumProb - dummy1
  # find the upper bound of the acceptance region
  SumProb2 <- 1
  yiNup <- 2*NN
  while(SumProb2 >= (1-(alpha/2))){
    dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNup)
    SumProb2 <- SumProb2 - dummy1
    yiNup <- yiNup - 1
  }
  yiNup <- yiNup + 1
  dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNup)
  SumProb2 <- SumProb2 + dummy1
  # Test if yiNobs is with in the accptance region
  test <- yiNobs >= yiNlow && yiNobs <= yiNup
  if(test){
    output <- paste(yiNobs, " IS within the accptance region", sep = "")
  } else {
    output <- paste(yiNobs, " IS NOT within the accptance region", sep = "")
  }
  out <- list(lowerbound = yiNlow, 
              upperbound = yiNup,
              result = output)
  return(out)
}

AcceptanceRegion(1000, 30, 0.625, 0.25, 50, 0.05)
