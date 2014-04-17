\name{EqAppr-class}
\Rdversion{1.1}
\docType{class}
\alias{EqAppr-class}
\alias{calfig,EqAppr,character-method}
\alias{getModnames,EqAppr-method}
\alias{getPruned,EqAppr-method}
\alias{getSens,EqAppr-method}
\alias{getUnpruned,EqAppr-method}
\alias{show,EqAppr-method}
\alias{show,ApprRes-method}

\title{Class \code{"EqAppr"}}
\description{
Manage the appraisal of an eQTL search
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EqAppr", ...)}, or via \code{buildEqAppr()}
}
\section{Slots}{
  \describe{
    \item{\code{meta}:}{Object of class \code{"ApprMeta"} basic descriptive information
about source analysis }
    \item{\code{sens}:}{Object of class \code{"ApprSens"} outputs of sensitivity analysis }
    \item{\code{pruned}:}{Object of class \code{"ApprRes"} outputs of appraisal for LD-pruned predictions }
    \item{\code{unpruned}:}{Object of class \code{"ApprRes"} outputs of general appraisal without LD pruning }
  }
}
\section{Methods}{
  \describe{
    \item{calfig}{\code{signature(x = "EqAppr", ind = "character")}: 
      Additional arguments can be specified:
       \describe{
          \item{ind}{index or name of model to be plotted}
          \item{hfudgetxt}{distance to move rendered fractions relative
        to bin x coordinate}
          \item{tickend}{
        maximum value at which axis tick mark will be plotted
        }
          \item{tickgap}{
        axis will have ticks at \code{seq(0, tickend, tickgap)}
        }
          \item{ylimin}{
        ylim setting for rendering
        }
          \item{xlimin}{
        xlim setting for rendering
        }
          \item{fraccex}{
        cex setting for fraction rendering
        }
          \item{fuselast}{
        if data are sparse in entries of high predicted probability,
        you can fuse the nearby cells up to the end -- pick fuselast=2
        for final 2 cells, 3 for final 3 and so on
        }
       }
     }
    \item{getModnames}{\code{signature(x = "EqAppr")}: obtain the list
of strings used to name different appraisal models }
    \item{getPruned}{\code{signature(x = "EqAppr")}: get the
ApprRes instance corresponding to the LD-pruned loci }
    \item{getUnpruned}{\code{signature(x = "EqAppr")}: get the
ApprRes instance for all loci in use }
    \item{getSens}{\code{signature(x = "EqAppr")}: get sensitivity
analysis results }
    \item{show}{\code{signature(object = "EqAppr")}: concise report }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("EqAppr")
}
\keyword{classes}
