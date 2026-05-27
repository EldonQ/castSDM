\name{cast_safe_ggsave}
\alias{cast_safe_ggsave}
\title{Safely Save a ggplot Object}
\usage{
cast_safe_ggsave(filename, plot = ggplot2::last_plot(), ...)
}
\arguments{
\item{filename}{Output filename.}

\item{plot}{Plot object.}

\item{...}{Additional arguments passed to \code{\link[ggplot2:ggsave]{ggplot2::ggsave()}}.}
}
\value{
Invisibly returns \code{filename}.
}
\description{
Wrapper around \code{\link[ggplot2:ggsave]{ggplot2::ggsave()}} that applies
castSDM's configured plotting font and chooses conservative graphics devices.
PNG output uses \pkg{ragg} when available and otherwise falls back to base PNG.
PDF output uses Cairo when available and otherwise base PDF with an explicit
font family.
}
\seealso{
\code{\link[=cast_set_plot_defaults]{cast_set_plot_defaults}}
}
