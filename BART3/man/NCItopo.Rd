\name{NCItopo}
\docType{data}
\alias{NCItopo}
\title{NCI AntiHIV Drug Screening}
\description{From the National Cancer Institute (NCI)
  Developmental Therapeutics Program.  The discovery
  of new drugs depends on developing compounds with minimal toxic effects.
  Statistical predictive toxicology relies on observed training data
  to learn the relationship between chemical structure features and
  toxicity response.  
}

\usage{data(NCItopo)}

\details{CEM cells are derived from T cells provided by a patient (with
  initals CEM) who suffered from acute lymphoblastic leukemia.
  \code{Potency} is the categorical response that measures whether a
  compound protects CEM cells from HIV-1 infection.  Presumably,
  compounds with \code{Potency=0} are inactive and \code{Potency>0} are
  active.  The structural features provided are the Topological
  Information Indices (TII).  TII are a common choice since they are
  easy to calculate; very sensitive to small changes in molecular
  structure; and do not depend on conformation of the molecule.  This is
  a data-frame with 29,374 observations and 262 variables.  Except for
  \code{Potency} and \code{Name}, there are 260 TII variables that are
  included. However, due to the vast number, TII covariates are not
  described further.  See the Todeschini and Consonni reference
  (pp.447-456) for more information.}

\references{
  Feng et al. (2003)
  Predictive toxicology: benchmarking molecular descriptors and
  statistical methods. \emph{J Chem Inf Comput Sci} \bold{43}, 1463-1470.
  
  Todeschini R, Consonni V. (2000)
  Handbook of Molecular Descriptors. Viley-VCH: Weinheim, Germany.
}

\keyword{datasets}
