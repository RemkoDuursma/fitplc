#' An example vulnerability curve
#' @description Percent loss conductivity as a function of water potential for three species.
#' @name stemvul
#' @docType data
#' @format
#' \describe{
#' \item{Species}{One of dpap, egran, ssay}
#' \item{Branch}{Replicate branch, multiple branches were measured for each species}
#' \item{MPa}{Xylem water potential (MPa)}
#' \item{PLC}{Percent loss conductivity}
#' \item{Cond}{Raw, unscaled conductivity of branch segment (units)}
#' }
NULL


# Conductivity
# dfr <- read.csv("test/stemvul-points.csv", sep=";")
# dfr$Cond <- 5.5 * (100 - dfr$PLC)/100
# dfr$Cond[dfr$Species == "dpap"] <- 7.8 * (100 - dfr$PLC[dfr$Species == "dpap"])/100
# dfr$Cond[dfr$Species == "egran"] <- 12.2 * (100 - dfr$PLC[dfr$Species == "egran"])/100
# dfr$Rep <- NULL
# dfr$Group <- NULL
# stemvul <- dfr
# save(stemvul, file="data/stemvul.rda")
