url <- "https://raw.githubusercontent.com/davidearn/London_smallpox/master/Data/London_smallpox.csv"
csv <- tempfile()
download.file(url, csv)

dd0 <- read.csv(csv)

dd1 <- transform(dd0,
                 from = as.Date(sprintf("%04d-%02d-%02d",
                                        year.from, month.from, day.from)),
                 to   = as.Date(sprintf("%04d-%02d-%02d",
                                        year.to  , month.to  , day.to  )))

dd2 <- dd1[rowSums(is.na(dd1[c("smpx", "acm", "birth")])) < 3L, ]

dd3 <- with(dd2,
            data.frame(from = .Date(as.integer(from)),
                       nday = .difftime(as.integer(to - from), "days"),
                       smallpox = smpx, allcauses = acm, births = birth,
                       row.names = NULL)) # guarantees automatic row names

## One discrepancy that we don't ignore:
table(dd3[["nday"]])
r <- 11366L
f3 <- .Date(as.integer(as.Date(c("1881-12-17", "1882-12-24", "1882-01-01"))))
with(dd3, {
stopifnot(exprs = {
	identical(sum(nday <= 0L), 1L)
	identical(which(nday <= 0L), r)
	identical(from[r + (-1L:1L)], f3)
})
})

dd4 <- dd3
dd4[r, "from"] <- dd4[r - 1L, "from"] + 7L
dd4[r, "nday"] <- 8L
table(dd4[["nday"]])

smallpox <- dd4
save(smallpox, file = "smallpox.rda", version = 3L, compress = "xz")
