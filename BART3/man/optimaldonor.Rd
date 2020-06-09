\name{optimaldonor}
\docType{data}
\alias{optimaldonor}
\title{Hematopoietic stem cell transplantation (HSCT) optimal donor}
\description{
  This data set contains the characteristics of hypothetical patients
  and their hypothetical donors.  With this data frame and the
  \code{predict} function, you can evaluate the potential improvement
  (or lack thereof) that a younger (older) donor would provide.
  See the demo \code{optimaldonor}.}
\usage{data(optimaldonor)}
\format{A data frame with 100 observations on the following 33 variables.
  \describe{
    \item{\code{agegp}}{Recipient age group (2=18-29, 3=30-39,
      4=40-49, 5=50-59, 6=60+)}
    \item{\code{sex}}{Recipient gender (1=Male, 2=Female)}
    \item{\code{cau}}{Caucasian (1=Yes, 0=No)}
    \item{\code{afa}}{African-American (1=Yes, 0=No)}
    \item{\code{api}}{Asian/Pacific Islander (1=Yes, 0=No)}
    \item{\code{nam}}{Native American (1=Yes, 0=No)}
    \item{\code{raceo}}{Other Race (1=Yes, 0=No)}
    \item{\code{raceu}}{Unknown Race (1=Yes, 0=No)}
    \item{\code{hisp}}{Hispanic (1=Yes, 0=No)}
    \item{\code{kpslt90}}{Karnofsky score <90 (1=Yes, 0=No)}
    \item{\code{kpsge90}}{Karnofsky score 90+ (1=Yes, 0=No)}
    \item{\code{kpsu}}{Karnofsky score unknown (1=Yes, 0=No)}
    \item{\code{aml}}{Acute myeloid leukemia at HSCT (1=Yes, 0=No)}
    \item{\code{all}}{Acute lymphoblastic leukemia at HSCT (1=Yes, 0=No)}
    \item{\code{cml}}{Chronic myeloid leukemia at HSCT (1=Yes, 0=No)}
    \item{\code{mds}}{Myelodysplastic syndrome at HSCT (1=Yes, 0=No)} 
    \item{\code{disstat}}{Disease status prior to HSCT (1=Early,
      2=Intermediate, 3=Advanced)}
    \item{\code{intdxtxgp}}{Time from diagnosis to HSCT (1= <6 months,
      2=6-12 months or missing, 3=12+ months)}
    \item{\code{intdxtxgpmi}}{Time from diagnosis to HSCT missing (1=Yes, 0=No)} 
    \item{\code{condint}}{Conditioning intensity (1=MAC, 2=RIC)}
    \item{\code{tbin}}{Total body irradiation (1=No, 0=Otherwise)}
    \item{\code{tbiy}}{Total body irradiation (1=Yes, 0=Otherwise)}
    \item{\code{tbiu}}{Total body irradiation (1=Missing, 0=Otherwise)}
    \item{\code{graftype}}{Graft type (1=Bone Marrow, 2=Peripheral Blood)}
    \item{\code{dnrage}}{Donor age in decades}  
    \item{\code{tac}}{Tacrolimus (FK506) based GVHD prophylaxis (1=Yes, 0=No)} 
    \item{\code{csa}}{Cyclosporine A (CSA) based GVHD prophylaxis (1=Yes, 0=No)}
    \item{\code{gvhdo}}{Non-FK506/Non-CSA based GVHD prophylaxis (1=Yes, 0=No)}
    \item{\code{invivo_tcd}}{In vivo T-cell depletion (1=Yes, 0=No)}
    \item{\code{yeargp}}{Year of HSCT (1=1999-2002, 2=2003-2006, 3=2007-2011,
      4=2012-2014, 5=2015-2016)}
    \item{\code{rcmvp}}{Recipient Cytomegalovirus Positive (1=Yes, 0=No)}
    \item{\code{rcmvn}}{Recipient Cytomegalovirus Negative (1=Yes, 0=No)}
    \item{\code{rcmvu}}{Recipient Cytomegalovirus Unknown (1=Yes, 0=No)}
  }
}
\details{Based on the collected data of hematopoietic stem cell
  transplantation outcomes derived from the registry compiled by
  the Center for International Blood and Marrow Research in
  association with the National Marrow Donor Program, we
  determined that, generally, the youngest available matched
  adult donor is preferable to the alternatives. We believe this
  is the current state-of-the-art.  However, as further
  data is collected, this evidence may need to be re-evaluated.
}

\keyword{datasets}
