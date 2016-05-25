`seemslog` <-
function(groups) {
    nologevidence = sum(sapply(groups, function(x) {
        ifelse(length(x) > 2, log(max(0, shapiro.test(x)$p.value)), 
            0)
    }))
    if (min(unlist(groups)) <= 0) {
        logevidence = -1e+200
    }
    else {
        logevidence = sum(sapply(groups, function(x) {
            ifelse(length(x) > 2, log(max(0, shapiro.test(log(x))$p.value)), 
                0)
        }))
    }
    (logevidence > nologevidence + 1)
}

