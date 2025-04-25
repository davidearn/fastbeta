## MJ: for reference only as very obfuscating to me ...

ss <- mgcv::s(times,
              bs = "ps",
              k = basis.dimension,
              m = c(basis.order - 2L, penalty.order))
str(ss)
sc <- mgcv::smooth.construct(object = ss,
                             data = data.frame(times),
                             knots = NULL)
str(sc)
sC <- mgcv::smoothCon(object = ss,
                      data = data.frame(times),
                      knots = NULL,
                      absorb.cons = FALSE,
                      scale.penalty = FALSE)[[1L]]
str(sC)
sR <- mgcv::smooth2random(object = sC,
                          vnames = "",
                          type = 2L)
str(sR)
