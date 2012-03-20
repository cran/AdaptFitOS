 ########## R function: aspmPenRead ##########

# For extracting information from a 
# penalized spline term.


"aspmPenReadOS" <-
  function (term) 
{
  arg.list <- substring(term, 3, (nchar(term) - 1))
  var.name <- break.string(arg.list, ",")[1]
  var.val <- eval(parse(text = var.name))    #   ,sys.frame(-3)
  out <- arg.searchOS(arg.list, "knots=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE) knots <- default.knots(var.val)
  if (present == TRUE)  knots <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "var.knots=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    var.knots <- NULL
  if (present == TRUE)     var.knots <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "adap=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    adap <- FALSE #TRUE
  if (present == TRUE)     adap <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "spar=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    spar <- NULL
  if (present == TRUE)     spar <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "adf=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    adf <- "miss"
  if (present == TRUE)     adf <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "basis=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    basis <- "os"
  if (present == TRUE)     basis <- spmArgRead(arg)$val
  out <- arg.searchOS(arg.list, "degree=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE){
    if (basis == "trunc.poly") degree <- 1
    else degree <- 3
  }
  if (present == TRUE) degree <- spmArgRead(arg)$val

  out <- arg.searchOS(arg.list, "var.basis=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    var.basis <- basis
  if (present == TRUE)     var.basis <- spmArgRead(arg)$val

  out <- arg.searchOS(arg.list, "var.degree=")
  present <- out$present
  arg <- out$arg
  if (present == FALSE)    var.degree <- degree
  if (present == TRUE)     var.degree <- spmArgRead(arg)$val

  return(list(name = var.name, var = var.val, knots = knots,
              adap = adap, var.knots = var.knots, var.basis=var.basis, var.degree=var.degree, spar = spar, adf = adf,
              degree = degree, basis = basis))
}
