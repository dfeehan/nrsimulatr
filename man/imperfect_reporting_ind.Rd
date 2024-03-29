% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporting_graph_ind.R
\name{imperfect_reporting_ind}
\alias{imperfect_reporting_ind}
\title{Create a set of parameters describing reporting behavior that is not perfectly accurate
and that varies by individual}
\usage{
imperfect_reporting_ind(params)
}
\arguments{
\item{params}{a list of params, which must contain tau}
}
\value{
a reporting parameter object
}
\description{
Create a set of parameters describing reporting behavior that is not perfectly accurate
and that varies by individual
}
\details{
For now, the only simulation parameter that makes a difference is
\eqn{\tau}.
}
