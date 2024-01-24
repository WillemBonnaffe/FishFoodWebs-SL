# Load the Rcpp package
library(Rcpp)

# Source the C++ code
sourceCpp("DEMCpp_v0.2.cpp")

# Define the bivariate Gaussian target distribution
dTarget <- function(theta) {
  mu <- c(0, 3)
  Sigma <- matrix(c(1, 0, 0, 1), nrow = 2)
  -0.5 * (theta - mu) %*% solve(Sigma) %*% (theta - mu)
}

# Set parameters
Theta_0 <- c(0, 0)
gamma <- 0.1
epsilon <- 0.1
nIt <- 900

# Generate primer chain
set.seed(123)
primerChain <- matrix(0, nrow = 100, ncol = length(Theta_0) + 1)
primerChain[, -1] <- matrix(rnorm(100 * length(Theta_0), mean = 0, sd = 1), ncol = length(Theta_0))
primerChain[, 1] <- apply(primerChain[, -1, drop = FALSE], 1, function(theta) dTarget(theta))
# primerChain = t(primerChain[1,])
head(primerChain)
# or use previous chain
primerChain = result$chainList

# Create argument list
argList <- list(
  dTarget = dTarget,
  # Theta_0 = Theta_0,
  # gamma = gamma,
  epsilon = epsilon,
  nIt = nIt,
  primerChain = primerChain
)

# Run the modified DEMCpp function
result <- DEMCpp(argList)

# Print the acceptance ratio and summary of the results
cat("Acceptance Ratio:", result$p, "\n")
summary(result$chainList)

## Check
tail(result$chainList[1:100,])
tail(primerChain)

## 
chain = result$chainList
par(mfrow=c(3,1))
plot(chain[,1],type="l")
plot(chain[,2],type="l")
plot(chain[,3],type="l")
par(mfrow=c(1,1))