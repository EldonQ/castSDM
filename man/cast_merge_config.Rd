\name{cast_merge_config}
\alias{cast_merge_config}
\title{Merge castSDM Configuration Lists}
\usage{
cast_merge_config(defaults, overrides = list())
}
\arguments{
\item{defaults}{Named list of default values.}

\item{overrides}{Named list of dataset-specific values.}
}
\value{
A named list.
}
\description{
Thin wrapper around \code{\link[utils:modifyList]{utils::modifyList()}}
with a clearer name for scripts.
}
\seealso{
\code{\link[=cast_default_config]{cast_default_config}},
\code{\link[=cast_run_from_config]{cast_run_from_config}}
}
