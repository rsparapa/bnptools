\name{bmx}
\alias{bmx}
\docType{data}
\title{NHANES 1999-2000 Body Measures and Demographics}

\description{ This data set was created from the National Health and
  Nutrition Examination Survey (NHANES) 1999-2000 Body Measures Exam and
  Demographics. To create growth charts, this data is restricted to
  3435 children aged 2 to 17.
  %% The relative uncertainty is derived
  %% from the CDC growth charts based on the width from the 2.5 to
  %% the 97.5 percentile; see \code{CDCweight} and \code{CDCheight}.
}

\usage{
data(bmx)
}

\format{
  \tabular{ll}{
    SEQN:\tab Sequence number\cr
    height:\tab BMXHT in cm\cr
    sex:\tab RIAGENDR 1=male, 2=female\cr
    age:\tab RIDAGEEX in years with fractions for months\cr
    weight:\tab BMXWT in kg\cr
    race:\tab RIDRETH2 1=White, 2=Black, 3=Hispanic\cr
    %%CDCwt:\tab Relative uncertainty for weight\cr
    %%CDCht:\tab Relative uncertainty for height\cr
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
