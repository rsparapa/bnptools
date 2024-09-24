\name{bmx}
\alias{bmx}
\docType{data}
\title{NHANES 1999-2000 Body Measures and Demographics}

\description{ This data set was created from the National Health and
  Nutrition Examination Survey (NHANES) 1999-2000 Body Measures Exam and
  Demographics. To create growth charts, this data is restricted to
  3435 children aged 2 to 17. The relative uncertainty is derived
  from the CDC growth charts based on the width from the 2.5 to
  the 97.5 percentile; see \code{CDCweight} and \code{CDCheight}.}

\usage{
data(bmx)
}

\format{
  \tabular{ll}{
    SEQN:\tab Sequence number\cr
    BMXHT:\tab Height in cm\cr
    RIAGENDR:\tab Gender: 1=male, 2=female\cr
    RIDAGEEX:\tab Age in years with fractions for months\cr
    BMXWT:\tab Weight in kg\cr
    RIDRETH2:\tab Race/ethnicity: 1=Non-Hispanic White,
    2=Non-Hispanic Black, 3=Hispanic\cr
    CDCwt:\tab Relative uncertainty for weight\cr
    CDCht:\tab Relative uncertainty for height\cr
  }
}

\references{
  National Health and Nutrition Examination Survey 
  1999-2000 Body Measures Exam.
  \url{https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/BMX.htm}

  National Health and Nutrition Examination Survey 
  1999-2000 Demographics.
  \url{https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/DEMO.htm}
}

\keyword{datasets}
