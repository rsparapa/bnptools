\name{hnscc}
\docType{data}
\alias{hnscc}
\title{Head and Neck Squamous Cell Carcinoma}

\description{

  Head and Neck Squamous Cell Carcinoma (HNSCC) is the sixth most common
	cancer worldwide.  Most HNSCC patients present with
	loco-regionally advanced disease, stage III/IV, requiring
	multimodal treatment including surgery, radiation, and/or
	chemotherapy. The standard of care for patients with resectable
	disease is salvage surgery. Risk factors for recurrence of HNSCC
	after salvage are extra-nodal extension (ENE), positive surgical
	margins, perineural invasion (PNI), lymphovascular invasion
	(LVI), and involvement of 2 or more lymph nodes. A single-arm
	phase II study was conducted for evaluating disease-free
	survival after 6 months of adjuvant nivolumab following salvage
	resection for locally recurrent HNSCC with the potential
	efficacy of PD-1 inhibitors (Leddon et al. 2022).  Historical
	controls extracted from electronic medical records at the same
	institution were used to determine treatment efficacy (Haque et
	al. 2019).

      }

      \usage{
      data(hnscc)

}

\format{A data frame with 105 subjects on the following 17 variables.
\describe{
\item{\code{DFS}}{an indicator for disease-free survival: 1=event, 0=censoring}
\item{\code{DFS_Time}}{time-to-event in months}
\item{\code{Age}}{age in years}
\item{\code{Gender}}{sex: a factor with levels \code{Female} and \code{Male}}
\item{\code{Race}}{race: a factor with levels \code{White}, \code{Black} and \code{Asian}}
\item{\code{Smoking}}{history of 10 or more pack-years: a factor with levels \code{Yes} and \code{No}}	
\item{\code{Alcohol}}{consumption of 5 or more drinks per week: a factor with levels \code{Yes} and \code{No}}	
\item{\code{Prior_Chemo}}{a factor with levels \code{Yes} and \code{No}}
\item{\code{Primary_Disease_Site}}{a factor with levels \code{Larynx}, \code{Oral Cavity} and \code{Oropharynx}}
\item{\code{ECS}}{extracapsular spread from lymph nodes: a factor with levels \code{Yes}, \code{No} and blank for missing}	
\item{\code{Positve_M}}{positive surgical margins: a factor with levels \code{Yes}, \code{No} and blank for missing}
\item{\code{PNI}}{perineural invasion: a factor with levels \code{Yes}, \code{No} and blank for missing}
\item{\code{LVI}}{lymphovascular invasion: a factor with levels \code{Yes}, \code{No} and blank for missing}
\item{\code{Risk_Stratification}}{a factor with levels \code{High Risk}
  (positive surgical margins and extranodal extension),
  \code{Not High Risk} and blank for missing}
\item{\code{Lymph_Nodes_Involved_2}}{involvement of >2 lymph nodes: a factor with levels \code{Yes}, \code{No} and blank for missing}
\item{\code{Study}}{\code{1} (trial population) and \code{0} (historical controls)}
\item{\code{diag_year}}{\code{1} and \code{2}}
}
}

\references{
  Haque S, Karivedu V, Riaz MK, Choi D, Roof L, Hassan SZ, Zhu Z,
  Jandarov R, Takiar V, Tang A, Wise-Draper T. (2019)
  "High-risk pathological features at the time of salvage surgery
  predict poor survival after definitive therapy in patients with
  head and neck squamous cell carcinoma".
  \emph{Oral Oncol}, \bold{88} pp. 9-15.
	
  Leddon JL, Gulati S, Haque S, Allen C, Palackdharry S, Mathews M, 
  Kurtzweil N, Riaz MK, Takiar V, Nagasaka M, Patil Y, Zender C, Tang A, 
  Cervenka B, McGrath J, Korn WM, Hinrichs BH, Jandarov R, Harun N, Sukari A, 
  Wise-Draper TM. (2022)
  "Phase II Trial of Adjuvant Nivolumab Following Salvage Resection 
  in Patients with Recurrent Squamous Cell Carcinoma of the
  Head and Neck". 
  \emph{Clin Cancer Res}, \bold{28(16)}, pp.=3464-3472. 
}

\examples{
  \dontrun{
      ## you can analyze this data set with BART as follows
      source(system.file('HNSCC/bart.R', package = 'BART3'))
      ## see the HNSCC sub-directory above for more examples 
    }
  }

\keyword{datasets}
