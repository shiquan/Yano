.normargPwm <- function(pwm, argname="pwm")
{
    if (!is.matrix(pwm) || !is.numeric(pwm))
        stop("'", argname, "' must be a numeric matrix")
    if (!identical(rownames(pwm), DNA_BASES))
        stop("'rownames(", argname, ")' must be the 4 DNA bases ('DNA_BASES')")
    if (!is.double(pwm))
        storage.mode(pwm) <- "double"
    if (any(is.na(pwm)))
        stop("'", argname, "' contains NAs")
    pwm
}

.normargMinScore <- function(min.score, pwm)
{
    if (!isSingleNumber(min.score) && !isSingleString(min.score))
        stop("'min.score' must be a single number or string")
    if (is.numeric(min.score)) {
        if (!is.double(min.score))
            storage.mode(min.score) <- "double"
        return(min.score)
    }
    nc <- nchar(min.score)
    if (substr(min.score, nc, nc) == "%")
        min.score <- substr(min.score, 1L, nc-1L)
    maxScore(pwm) * as.double(min.score) / 100.00
}

setGeneric("matchMOTIF", signature="subject",
    function(pwm, subject, min.score="80%", ...)
        standardGeneric("matchMOTIF")
    )

#' @import Biostrings
#' 
setMethod("matchMOTIF", "DNAString",
          function(pwm, subject, min.score="80%") {
            ## checking 'pwm'
            pwm <- .normargPwm(pwm)
            ## checking 'min.score'
            min.score <- .normargMinScore(min.score, pwm)
            base_codes <- xscodes(subject, baseOnly=TRUE)
            C_ans <- .Call2("match_PWM_fast", pwm, subject, min.score, base_codes, PACKAGE = "Yano")
            return(C_ans)
          }
          )


setGeneric("RBP_matchmotifs",
           function(pwms, subject, genome, min.score, ...)
             standardGeneric("RBP_matchmotifs")
           )

setMethod("RBP_matchmotifs", signature(pwms="PWMatrixList", subject="GenomicRanges", genome="BSgenome"),
          function(pwms,
                   subject,
                   genome,
                   min.score = "80%"
                   ) {

            
          }
          )
