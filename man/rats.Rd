\name{rats}
\alias{rats}
\title{Growth curves}
\description{
  The \code{rats} data frame has 135 rows and 4 columns of the
  change in weight measured over time for rats assigned to different treatment groups.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{
       grouping variable.
    }
    \item{time}{
      time (week) of measurement (0, 1, 2, 3, 4).
    }
    \item{trt}{
      treatment group (1, 2, 3).
    }
    \item{y}{
      weight (grams)
    }
  }
}
\details{
  In a weight gain experiment, 30 rats were randomly assigned to three treatment groups: treatment 1, a control (no additive); treatments 2 and 3, which consisted of two different additives (thiouracil and thyroxin respectively) to the rats drinking water (Box, 1950). Weight (grams) of the rats was measured at baseline (week 0) and at weeks 1, 2, 3, and 4. Data on three of the 10 rats from the thyroxin group were subsequently removed due to an accident at the beginning of the study.
}
\source{
  G. E. P. Box, Problems in the analysis of growth and wear curves, Biometrics 6 (4) (1950) 362-389.
}

\keyword{datasets}
