\name{nlsym}
\docType{data}
\alias{nlsym}
\title{Monetary returns of post-secondary education}

\description{This is the data set from the well-known study by David
  Card.  This data frame contains 3613 observations on men in 1976: a
  cross-section of the US National Longitudinal Survey of Young Men
  (NLSYM).  NLSYM began in 1966 with 5525 men aged 14:24 and continued
  with follow-up surveys through 1981.}

\usage{data(nlsym)}

\format{A data frame with 3613 observations on the following 53 variables.
  \describe{
    \item{\code{id}}{sequential identifier from 1 to 5225}
    \item{\code{nearc2}}{grew up near 2-year college}
    \item{\code{nearc4}}{grew up near 4-year college}
    \item{\code{nearc4a}}{grew up near 4-year public college}
    \item{\code{nearc4b}}{grew up near 4-year priv college}
    \item{\code{ed76}}{years of education in 1976}
    \item{\code{ed66}}{years of education in 1966}
    \item{\code{age76}}{age in 1976}
    \item{\code{daded}}{education of father: missing data imputed by mean}
    \item{\code{nodaded}}{\code{daded} imputed}
    \item{\code{momed}}{education of mother: missing data imputed by mean}
    \item{\code{nomomed}}{\code{momed} imputed}
    \item{\code{weight}}{NLS weight for 1976 cross-section}
    \item{\code{momdad14}}{lived with mom and dad at age 14}
    \item{\code{sinmom14}}{lived with single mom age 14}
    \item{\code{step14}}{lived with step-parent age 14}
    \item{\code{reg661}}{region=1 in 1966: Census Division 1, New England}
    \item{\code{reg662}}{region=2 in 1966: Census Division 2, Middle Atlantic}
    \item{\code{reg663}}{region=3 in 1966: Census Division 3, East North Central}
    \item{\code{reg664}}{region=4 in 1966: Census Division 4, West North Central}
    \item{\code{reg665}}{region=5 in 1966: Census Division 5, South Atlantic}
    \item{\code{reg666}}{region=6 in 1966: Census Division 6, East South Central}
    \item{\code{reg667}}{region=7 in 1966: Census Division 7, West South Central}
    \item{\code{reg668}}{region=8 in 1966: Census Division 8, Mountain}
    \item{\code{reg669}}{region=9 in 1966: Census Division 9, Pacific}
    \item{\code{south66}}{lived in south in 1966: Census Divisions 5, 6 and 7}
    \item{\code{work76}}{worked in 1976}
    \item{\code{work78}}{worked in 1978}
    \item{\code{lwage76}}{log wage in 1976: outliers trimmed}
    \item{\code{lwage78}}{log wage in 1978: outliers trimmed}
    \item{\code{famed}}{mom-dad education class 1-9}
    \item{\code{black}}{African-American}
    \item{\code{smsa76r}}{SMSA in 1976}
    \item{\code{smsa78r}}{SMSA in 1978}
    \item{\code{reg76r}}{in south in 1976}
    \item{\code{reg78r}}{in south in 1978}
    \item{\code{reg80r}}{in south in 1980}
    \item{\code{smsa66r}}{in smsa in 1966}
    \item{\code{wage76}}{raw wage: cents per hour 1976}
    \item{\code{wage78}}{raw wage: cents per hour 1978}
    \item{\code{wage80}}{raw wage: cents per hour 1980}
    \item{\code{noint78}}{no interview in 1978}
    \item{\code{noint80}}{no interview in 1980}
    \item{\code{enroll76}}{enrolled in 1976}
    \item{\code{enroll78}}{enrolled in 1978}
    \item{\code{enroll80}}{enrolled in 1980}
    \item{\code{kww}}{1966 interview included a 28-item test called the "Knowledge of the World of Work"}
    \item{\code{iq}}{a normed IQ score}
    \item{\code{marsta76}}{marital status in 1976: 1=married}
    \item{\code{marsta78}}{marital status in 1976: 1=married}
    \item{\code{marsta80}}{marital status in 1976: 1=married}
    \item{\code{libcrd14}}{library card in home at age 14}
    \item{\code{exp76}}{experience in 1976}
  }
}

\references{
  Card, D. (1993) Using geographic variation in college proximity to
  estimate the return to schooling.
  \emph{National Bureau of Economic Research}. \bold{w4483}.
  \url{https://www.nber.org/papers/w4483}
}
