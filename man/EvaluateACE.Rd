% Generated by roxygen2 (4.0.2): do not edit by hand
\name{EvaluateACE}
\alias{EvaluateACE}
\title{MIREX ACE Evaluation}
\usage{
EvaluateACE(results.directory, analysis = c("SeventhsBass", "Root", "MajMin",
  "MajMinBass", "Sevenths", "Segmentation"), old.style = FALSE,
  adjust = "fdr")
}
\arguments{
\item{results.directory}{directory containing the results in Johan
Pauwels's output format}

\item{analysis}{chord vocabulary or segmentation to load}

\item{old.style}{whether to use old-style (Friedman) evaluation}

\item{adjust}{correction method for multiple comparisons (see \code{p.adjust})}
}
\description{
\code{EvaluateACE} runs a comparative evaluation of MIREX audio
chord estimation results for a particular chord vocabulary or
segmentation from the output of Johan Pauwels's evaluation scripts.
}

