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
if (FALSE) {
## For a smaller image:
attr(dd3, "nday") <- local({
    nday <- dd3[["nday"]]
    i <- which(nday != 7L)
    cbind(i, nday[i], deparse.level = 0L)
})
dd3[["nday"]] <- NULL
}

smallpox <- dd3
save(smallpox, file = "smallpox.rda", version = 3L, compress = "xz")
