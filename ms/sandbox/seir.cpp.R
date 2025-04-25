source("seir.R"); stopifnot(is.list(data))
dll <- "seir"
cpp <- paste0(dll, ".cpp")
lib <- paste0(dll, .Platform[["dynlib.ext"]])
TMB::compile(cpp)
dyn.load(lib)
parameters <- list(log_sd = 0,
                   log_size = 0,
                   b0 = double(data$R0),
                   b1 = double(data$R1))
obj <- TMB::MakeADFun(data = data, parameters = parameters, DLL = dll)
