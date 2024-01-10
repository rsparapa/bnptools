\name{xdm20.train}
\alias{xdm20.train}
\alias{xdm.train}
\docType{data}
\title{
A real data example for \code{recur.bart}. 
}
\description{
   
  A matrix containing a 20\% random sample as well as the whole 100\% of
  the training set for a real data example of recurrent events survival
  analysis.  There are 100 patients in the cohort: 50 in the training
  set and 50 in the testing set.  %The full data set, \code{xdm.train},
  %can be obtained online at
%\url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.tar}.
%\url{https://www.mcw.edu/-/media/MCW/Departments/Biostatistics/Tech-Report-64-April-2017-data-file.tar?la=en}.
There are 488 patients in the full cohort: 235 in the training set and
253 in the testing set.  See the Reference below (and the References
therein) for more detailed information; a brief synopsis follows.

We explored the hospital admissions for a cohort of patients with
diabetes cared for by the Froedtert and Medical College of Wisconsin
health network.  These patients were identified via their Electronic
Health Records (EHR) which include vital signs, diagnoses, procedures,
laboratory values, pharmacy orders and billing data.  This human
subjects research and de-identified data release was approved by the
Medical College of Wisconsin and Froedtert Hospital joint
Institutional Review Board.  To maintain patient privacy, roughly one
fourth of patients were randomly sampled for inclusion as well as other
de-identification procedures.

We identified likely incident diabetes mellitus type 2 patients by
tabulating their first diagnosis code of primary diabetes (ICD-9 codes
250.x0 and 250.x2) in 2006 or 2007, i.e., no such codes were found for
these patients prior to 2006 for as far back as each their records
go which is variable.  We restricted the population to adults aged 21 to
90 by 01/01/2008.  Among the patients treated in this health system, the
vast majority were racially self-identified as either white or black so
our inclusion criteria is restricted to these groups.  Since our
interest is in patients with primary diabetes, we excluded those
patients who were diagnosed with either secondary diabetes or
gestational diabetes.

For this cohort, we identified every hospital admission between
01/01/2008 and 12/31/2012.  For convenience, follow-up begins on
01/01/2008, rather than from each their actual incident diagnosis
date which varied over the preceding 2 years.  Following all patients
concurrently allows us to temporally adapt, via our model, for
seasonal/epidemic hospital admissions such as the H1N1 influenza
outbreak in the US from April to June 2009.

We investigated the following risk factors: gender, race, age, insurance
status (commercial, government or other), diabetes therapy (insulin,
metformin and/or sulfonylurea), health care charges, relative value
units (RVU), vital signs, laboratory values, comorbidity/complication
diagnoses and procedures/surgeries (we will refer to vital signs and
laboratory values collectively as signs; and comorbidity/complication
diagnoses and procedures/surgeries collectively as conditions).  In
total, we considered 85 covariates of which 82 are external factors as
described above and three are temporal factors: time, \eqn{t}, the
counting process, \eqn{N_i(t-)}, and the sojourn time, \eqn{v_i(t)}.
Among these potential predictors only gender, race and age are
time-independent.  The rest are defined as last value carried forward.

For insulin, metformin and sulfonylurea, we only had access to
prescription orders (rather than prescription fills) and self-reported
current status of prescription therapy during clinic office visits.
Since, generally, orders are only required after every three fills, and
each fill can be for up to 90 days, we define insulin, metformin and
sulfonylurea as binary indicators which are one if there exists an order
or current status indication within the prior 270 days; otherwise zero.

Health care charges and relative value units (RVU) are measures related
to the services and procedures delivered.  However, they are so closely
related that recent charges/RVUs are of no practical value in this
analysis.  For example, just prior to a hospital admission on
a non-emergent basis, they often have a series of diagnostic tests and
imaging.  Similarly, for an emergent admission, the patient is often
seen in the emergency department just prior to admission where similar
services are conducted.  We do not consider these charges/RVUs
predictive of an admission because we are interested in identifying
preventive opportunities.  Therefore, we investigate charges/RVUs that
are the sum total of the following moving windows of days prior to any
given date: 31 to 90, 91 to 180, 181 to 300.

For many patients, some signs were not available for a given date so
they were imputed; similarly, if a sign was not observed within the last
180 days, then it was imputed (except for height which never expires,
weight extended to 365 days and body mass index which is a deterministic
function of the two).  We utilized the Sequential BART missing
imputation method.  However, instead of creating several imputed data
sets, we imputed a new sign at each date when it was missing, i.e., in
order to properly address uncertainty within one data set, a new value
was imputed for each date that it was missing and never carried forward.

Conditions are binary indicators which are zero until the date of the
first coding and then they are one from then on.  Based on clinical
rationale, we identified 26 conditions (23 comorbidities and 3
procedures/surgeries) which are potential risk factors for a hospital
admission many of which are possible complications of diabetes; besides
clinical merit, these conditions were chosen since they are present in
more than just a few subjects so that they may be informative.
Similarly, we employed 15 general conditions which are the Charlson
diagnoses and 18 general conditions from the RxRisk adult diagnoses
which are defined by prescription orders.  Seven conditions are a
composite of diagnosis codes and prescription orders.
}
\format{A matrix with 84 columns of, mainly, binary indicators.
\describe{
  \item{\code{t}}{time in days}
  \item{\code{v}}{sojourn time since last event in days}
  \item{\code{n}}{number of previous events}
  \item{\code{age}}{age in years}
  \item{\code{sex_c}}{1:F, 2:M}
  \item{\code{race1}}{1:W, 2:B}
  \item{\code{insur}}{Insurance status 0:Government aged 65+, 1:Government aged<65, 2:Commercial aged<65, 3:Other aged<65}
\item{\code{ins270}}{Insulin prescription within last 270 days}
\item{\code{met270}}{Metformin prescription within last 270 days}
\item{\code{sulf270}}{Sulfonylurea prescription within last 270 days}
\item{\code{charges31_90}}{Health care charges 31 to 90 days prior}
\item{\code{charges91_180}}{Health care charges 91 to 180 days prior}
\item{\code{charges181_300}}{Health care charges 181 to 300 days prior}
\item{\code{rvus31_90}}{Relative value units 31 to 90 days prior}
\item{\code{rvus91_180}}{Relative value units 91 to 180 days prior}
\item{\code{rvus181_300}}{Relative value units 181 to 300 days prior}
\item{\code{keto}}{Ketoacidosis}
\item{\code{hypergly}}{Hyperglycemia}
\item{\code{hyperins}}{Hyperinsulinism}
\item{\code{ckd}}{Chronic kidney disease}
\item{\code{nephropath}}{Nephropathy}
\item{\code{kidfail}}{Kidney failure}
\item{\code{dialy}}{Dialysis}
\item{\code{kidtrans}}{Kidney transplant}
\item{\code{retinopath}}{Retinopathy}
\item{\code{blind}}{Blindness}
\item{\code{neuropath}}{Neuropathy}
\item{\code{dfoot}}{Diabetic foot}
\item{\code{enceph}}{Encephalopathy}
\item{\code{gang}}{Gangrene}
\item{\code{cmyopath}}{Cardiomyopathy}
\item{\code{cad}}{Coronary artery disease}
\item{\code{hypertsn}}{Hypertension}
\item{\code{afib}}{Atrial fibrillation}
\item{\code{ulcer}}{Diabetic ulceration}
\item{\code{dmbone}}{Diabetic bone changes}
\item{\code{apnea}}{Sleep apnea}
\item{\code{depress}}{Depression}
\item{\code{cci_mi}}{Myocardial infarction}
\item{\code{cci_chf}}{Congestive heart failure}
\item{\code{cci_pvd}}{Peripheral vascular disease}
\item{\code{cci_stroke}}{Cerebrovascular disease}
\item{\code{cci_dement}}{Dementia}
\item{\code{cci_cpd}}{Chronic pulmonary disease}
\item{\code{cci_rhd}}{Rheumatic disease}
\item{\code{cci_pud}}{Peptic ulcer disease}
\item{\code{cci_ldmild}}{Mild liver disease}
\item{\code{cci_ldsev}}{Moderate/severe liver disease}
\item{\code{cci_dmw}}{Diabetes with chronic complications}
\item{\code{cci_plegia}}{Hemiplegia/paraplegia}
\item{\code{cci_rend}}{Renal disease}
\item{\code{cci_malig}}{Malignancy}
\item{\code{cci_tumor}}{Metastatic solid tumor}
\item{\code{barsurg}}{Bariatric surgery}
\item{\code{cabg}}{Coronary artery bypass graft}
\item{\code{ptca}}{Percutaneous transluminal coronary angioplasty}
\item{\code{nutr}}{Medical nutrition therapy}
\item{\code{asthma}}{Asthma}
\item{\code{hyperlip}}{Hyperlipidemia}
\item{\code{anxiety}}{Anxiety and tension}
\item{\code{psycho}}{Psychosis}
\item{\code{epilepsy}}{Epilepsy}
\item{\code{hd}}{Heart disease}
\item{\code{gout}}{Gout}
\item{\code{cd}}{Cardiac disease}
\item{\code{thyroid}}{Thyroid disorder}
\item{\code{tb}}{Tuberculosis}
\item{\code{pd}}{Parkinson''s disease}
\item{\code{bp_systolic}}{Systolic blood pressure}
\item{\code{bp_diastolic}}{Diastolic blood pressure}
\item{\code{bmi}}{Body mass index}
\item{\code{weight}}{Body weight}
\item{\code{height}}{Body height}
\item{\code{temperature}}{Body temperature}
\item{\code{pulse}}{Heart rate}
\item{\code{respirations}}{Respiratory rate}
\item{\code{a1c}}{Glycohemoglobin A1c}
\item{\code{crea}}{Creatinine}
\item{\code{sod}}{Sodium}
\item{\code{pot}}{Potassium}
\item{\code{bcar}}{Bicarbonate}
\item{\code{bun}}{Blood urea nitrogen}
\item{\code{ca}}{Calcium}
\item{\code{cl}}{Chloride}
}
}
\usage{
data(xdm20.train)
data(xdm.train)
}
\seealso{xdm20.test}
\references{
Sparapani, Rein, Tarima, Jackson, Meurer (2020).  Non-parametric recurrent
events analysis with BART and an application to the hospital admissions
of patients with diabetes. \emph{Biostatistics} doi:10.1093/biostatistics/kxy032
%% Sparapani R., Rein L., Tarima S., Jackson T., Meurer J. (2017)
%% Nonparametric recurrent events analysis with BART
%% and an application to the hospital admissions of patients with diabetes.
%% \emph{MCW Biostatistics Technical Report 64}
%% \url{http://www.mcw.edu/FileLibrary/Groups/Biostatistics/TechReports/TechReports5175/tr064.pdf}.

}
\examples{
data(xdm20.train)
head(xdm20.train[ , 1:10])
}
\keyword{datasets}
