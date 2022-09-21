#' linreg
#' 
#' A Reference Class that creates a linear regression based on a given
#' formula with p independent variables and a given data set with n 
#' observations. Allows users to view plots and data summaries based on this 
#' regression.
#' 
#' @field regdata The name of the data set being used. A character string.
#' @field regformula The formula being used for the regression. A character
#' string.
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
#' @return An object of class "linreg".
#' 
#' @export linreg
#' @exportClass linreg
#' @importFrom methods new setRefClass 
#' @import ggplot2


linreg <- setRefClass("linreg", 
                      fields = list(regdata = "character",
                                    regformula = "character",
                                    coefficients = "matrix", 
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
    
    "Creates the object, executes the regression calculations based on formula
    and data, and saves the values under the respective fields of the object."
    
    # First, we must establish our model matrix, our dependent variable, values
    # for n and p, and our y matrix to do our calculations:
    
    .self$regdata <- deparse(substitute(data))  # our data name
    .self$regformula <- format(formula)         # our formula name
    xmatrix <- model.matrix(formula, data)      # our model matrix
    y <- all.vars(formula)[1]                   # y, our dependent variable
    n <- nrow(data)                             # n, or number of observations
    p <- length(all.vars(formula))              # p, or number of parameters
    ymatrix <- matrix(data[, y])                # matrix with the values of y
    
    for (i in 2:length(all.vars(formula))){
      string = all.vars(formula)[i]
      
      if (is(data[1, string], "factor")){
        adds = length(unique(data[, string])) - 2
        p = p + adds
      }
    }
    
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
    for(i in 1:p){
      bhatvarvector <- append(bhatvarvector, bhatvarmatrix[i,i])
    }
    
    .self$coeffvariance <- matrix(bhatvarvector, 
                                  nrow = p)
    rownames(coeffvariance) <<- rownames(bhatvarmatrix)
    
    # This format allows us to easily create a matrix of the standard errors
    # and use that to get our t-values, and by extension our p-values:
    
    bhatsterr <- sqrt(coeffvariance)
    .self$tvalues <- coefficients / bhatsterr
    
    .self$pvalues <- pt(abs(tvalues), degrees, lower.tail = FALSE)
    
    # Finally, we place our summary statistics in a matrix so that we can
    # easily print it in the summary() method:
    
    .self$sumstats <- round(cbind(coefficients, sqrt(coeffvariance), 
                                  tvalues, pvalues), 2)
    
    colnames(sumstats) <<- c("Estimate", "Std. Error", "T Value", "P Value")
  },
  
  
  print = function(){
    
    "Prints the coefficient estimates for the regression."
    
    cat("linreg(formula = ", regformula, ", ", "data = ", regdata, ")",
        sep = "")
    cat("\n", "\n", rownames(coefficients))
    cat("\n", as.vector(coefficients))
  },
  
  plot = function(){
    
    "Shows 2 plots based on the regression: Residuals vs Fitted Values and 
    Scale - Location"
    
    theme_liu <- function(){
      
      ggplot2::theme_minimal() %+replace%
        
        ggplot2::theme(panel.border = element_rect(colour = "Light Grey",
                                                   fill = NA),
                       plot.margin = margin(30, 30, 20, 20),
                       plot.title = element_text(size = 18,
                                                 margin = margin(0,0,5,0),
                                                 vjust = 3,
                                                 hjust = 0),
                       plot.subtitle = element_text(size = 12,
                                                    vjust = 4,
                                                    margin = margin(5,0,5,0),
                                                    hjust = 0),
                       axis.title.x = element_text(size = 12),
                       axis.title.y = element_text(size = 12,
                                                   angle = 90),
                       axis.text.x = element_text(size = 10,
                                                  margin = margin(5,0,10,0),
                                                  colour = "black"),
                       axis.text.y = element_text(size = 10,
                                                  hjust = 1,
                                                  margin = margin(0,5,0,10),
                                                  colour = "black"),
                       axis.ticks = element_line(colour = "light grey"),
                       axis.ticks.length = unit(10, "pt"))
    }
    
    median <- median(residuals)
    g_plot_rvf <- 
      ggplot2::ggplot(mapping = aes(x = fittedvalues, y = residuals)) + 
      ggplot2::geom_point(size = 2, colour = "red") + 
      ggplot2::stat_summary(fun = "median", colour = "red", geom = "line") + 
      ggplot2::scale_x_continuous(limits = c(min(fittedvalues - 1),
                                             max(fittedvalues + 1))) + 
      ggplot2::scale_y_continuous(limits = c(min(residuals - sd(residuals)),
                                             max(residuals + sd(residuals)))) + 
      ggplot2::labs(title = "Residuals vs. Fitted Values", 
                    subtitle = paste("Regression:", regformula), 
                    x = "Fitted Values", 
                    y = "Residuals") + 
      theme_liu()
      
    
    s_res <- sqrt(abs(residuals/(sqrt(residvariance))))
    g_plot_srvf <- 
      ggplot2::ggplot(mapping = aes(x = fittedvalues, y = s_res)) +  
      ggplot2::geom_point(size = 2, colour = "red") + 
      ggplot2::stat_summary(fun = "mean", colour = "red", geom = "line") + 
      ggplot2::scale_x_continuous(limits = c(min(fittedvalues - sd(s_res)),
                                             max(fittedvalues + sd(s_res)))) + 
      ggplot2::scale_y_continuous(limits = c(min(s_res - .2),
                                             max(s_res + .2))) + 
      ggplot2::labs(title = "Scale - Location", 
                    subtitle = paste("Regression:", regformula), 
                    x = "Fitted Values", 
                    y = expression(sqrt(abs("Standardized Residuals")))) + 
      theme_liu()
    
    plots <- list(g_plot_rvf, g_plot_srvf)
    plots
    
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
    
    for(i in 1:nrow(sumstats)){
      stars = NA
      if (sumstats[i,4] < 0.001){
        stars = "***"
      } else if (sumstats[i,4] < 0.01){
        stars = "**"
      } else if (sumstats[i,4] < 0.05){
        stars = "*"
      } else if (sumstats[i,4] < 0.1){
        stars = "."
      } else {
        stars = " "
      }
      
      cat(rownames(sumstats)[i], sumstats[i,1], sumstats[i,2], sumstats[i,3],
          sumstats[i,4], stars, "\n", sep = " ")
    }
    
    cat("---", "\n", "Residual standard error:", round(sqrt(residvariance), 2), 
        "on", degrees, "degrees of freedom")
  }
)

