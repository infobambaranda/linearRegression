#' linreg
#' 
#' A Reference Class that creates a linear regression based on a given
#' formula with p independent variables and a given data set with n 
#' observations. Allows users to view plots and data summaries based on this 
#' regression.
#' 
#' @field coefficients Coefficient estimates for each independent variable. A
#' matrix with one column and p+1 rows.
#' @field fittedvalues Fitted values of dependent variable for each observation.
#' A matrix with one column and n rows.
#' @field residuals Residual values for dependent variable for each observation.
#' A matrix with one column and n rows.
#' @field degrees Degrees of Freedom. A length-one numeric vector.
#' @field residvariance Residual Variance. A length-one numeric vector.
#' @field coeffvariance Coefficient Variance for each independent variable. A 
#' matrix with one column and p+1 rows.
#' @field tvalues T-values for each coefficient estimate. A matrix with one 
#' column and p+1 rows.
#' @field pvalues P-values for each coefficient estimate. A matrix with one 
#' column and p+1 rows.
#' @field sumstats Summary statistics for the regression (coefficient estimates,
#' standard errors, t-values, and p-values). A matrix with with four columns
#' and p+1 rows.
#' 

linreg <- setRefClass("linreg", 
                      fields = list(coefficients = "matrix", 
                                    fittedvalues = "matrix", 
                                    residuals = "matrix", 
                                    degrees = "numeric", 
                                    residvariance = "numeric",
                                    coeffvariance = "matrix",
                                    tvalues = "matrix",
                                    pvalues = "matrix",
                                    sumstats = "matrix"))

linreg$methods(
  
  initialize = function(formula, data){
    
    "Creates the object, executes the regression calculations, and saves the 
    values under the respective fields."
    
    # First, we must establish our model matrix, our dependent variable, values
    # for n and p, and our y matrix to do our calculations:
    
    xmatrix <- model.matrix(formula, data)  # our model matrix
    y <- all.vars(formula)[1]               # y, our dependent variable
    n <- nrow(data)                         # n, or number of observations
    p <- length(all.vars(formula))          # p, or number of parameters
    ymatrix <- matrix(data[, y])            # matrix with the values of y
    
    
    # Next, we use these values to calculate the regression coefficients,
    # fitted values, residuals, degrees of freedom, residual variance,
    # coefficient variance, standard error, t-values, and p-values of 
    # the regression:
    
    # Regression coefficients
    .self$coefficients <- 
      solve((t(xmatrix) %*% xmatrix)) %*% t(xmatrix) %*% ymatrix
    
    # Fitted values
    .self$fittedvalues <- xmatrix %*% coefficients
    
    # Residuals matrix
    .self$residuals <- ymatrix - fittedvalues
    
    # Degrees of Freedom
    .self$degrees <- n - p
    
    # Residual variance
    .self$residvariance <- as.numeric((t(residuals) %*% residuals) / degrees)
    
    # To get our coefficient variances in a format we can use, we take the
    # values in the diagonal of the matrix produced below and put these values
    # in a new matrix with one column:
    
    bhatvarmatrix <- residvariance * solve((t(xmatrix) %*% xmatrix))
    
    bhatvarvector <- c()
    for(i in 1:length(all.vars(formula))){
      bhatvarvector <- append(bhatvarvector, bhatvarmatrix[i,i])
    }
    
    .self$coeffvariance <- matrix(bhatvarvector, 
                                  nrow = length(all.vars(formula)))
    rownames(coeffvariance) <<- rownames(bhatvarmatrix)
    
    # This format allows us to easily create a matrix of the standard errors
    # and use that to get our t-values, and by extension our p-values:
    
    bhatsterr <- sqrt(coeffvariance)
    .self$tvalues <- coefficients / bhatsterr
    
    .self$pvalues <- pt(abs(tvalues), degrees, lower.tail = FALSE)
    
    # Finally, we place our summary statistics in a matrix so that we can
    # easily print it in the summary() method:
    
    .self$sumstats <- cbind(coefficients, sqrt(coeffvariance), tvalues, pvalues)
    colnames(sumstats) <<- c("Estimate", "Std. Error", "T Value", "P Value")
  },
  
  
  print = function(){
    
    "Prints the coefficient estimates for the regression."
    
    cat("Coefficients:", "\n")
    t(testreg$coefficients)[1,]
  },
  
  
  resid = function(){
    
    "Shows the residual values for each observation in the data set."
    
    residuals
  },
  
  pred = function(){
    
    "Shows the fitted/predicted values for each observation in the data set."
    
    fittedvalues
  },
  
  coef = function(){
    
    "Shows the coefficient estimates as a named vector."
    
    coefvector <- as.vector(coefficients)
    names(coefvector) <- rownames(coefficients)
    coefvector
  },
  
  summary = function(){
    
    "Shows the summary statistics for the regression, including the coefficient
    estimates, standard errors, t-values, p-values, residual standard error, 
    and degrees of freedom."
    
    cat("Coefficients:", "\n")
    methods::show(sumstats)
    cat("---", "\n", "Residual Standard Error:", sqrt(residvariance), "on",
        degrees, "degrees of freedom")
  }
)

