\name{cast_set_plot_defaults}
\alias{cast_set_plot_defaults}
\title{Configure Safe Plot Fonts}
\usage{
cast_set_plot_defaults(family = "Arial")
}
\arguments{
\item{family}{Character. Font family passed to ggplot2. Default \code{"Arial"}.}
}
\value{
Invisibly returns \code{family}.
}
\description{
Sets the plotting font used by castSDM figures. The package default is
\code{"Arial"}; on Windows the family is registered with
\code{\link[grDevices:windowsFonts]{grDevices::windowsFont()}} when available.
Use \code{\link[=cast_safe_ggsave]{cast_safe_ggsave}} for export so PNG/PDF
devices handle fonts consistently.
}
\seealso{
\code{\link[=cast_safe_ggsave]{cast_safe_ggsave}}
}
