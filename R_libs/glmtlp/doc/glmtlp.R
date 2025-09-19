## ----include=FALSE------------------------------------------------------------
# to control the output
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
  lines <- options$output.lines
  if (is.null(lines)) {
    return(hook_output(x, options))  # pass to default hook
  }
  x <- unlist(strsplit(x, "\n"))
  more <- "..."
  if (length(lines)==1) {        # first n lines
    if (length(x) > lines) {
      # truncate the output, but add ....
      x <- c(head(x, lines), more)
    }
  } else {
    x <- c(more, x[lines], more)
  }
  # paste these lines together
  x <- paste(c(x, ""), collapse = "\n")
  hook_output(x, options)
})

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("glmtlp")

## -----------------------------------------------------------------------------
library(glmtlp)

## -----------------------------------------------------------------------------
data(gau_data)
X <- gau_data$X
y <- gau_data$y

## -----------------------------------------------------------------------------
fit <- glmtlp(X, y, family = "gaussian", penalty = "tlp")
fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
fit3 <- glmtlp(X, y, family = "gaussian", penalty = "l1")

## -----------------------------------------------------------------------------
plot(fit, xvar = "lambda")

## ----output.lines = 1:10------------------------------------------------------
coef(fit)
coef(fit, lambda = 0.1)

## -----------------------------------------------------------------------------
predict(fit, X[1:5, ], lambda = 0.1)
predict(fit, X[1:5, ], which = 10) # the 10th lambda in the lambda sequence

## -----------------------------------------------------------------------------
cv.fit <- cv.glmtlp(X, y, family = "gaussian", penalty = "tlp")

## -----------------------------------------------------------------------------
plot(cv.fit)

## ----output.lines = 10--------------------------------------------------------
coef(cv.fit)
predict(cv.fit, X[1:5, ])

