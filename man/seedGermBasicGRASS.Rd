\name{seedGermBasicGRASS}
\alias{seedGermBasicGRASS}
\title{Basic seed germination}
\usage{
  seedGermBasicGRASS(input, disturbance, suitability,
    pGermination, pGerminationDep, seedlToDcCuts,
    dcClasses, pGerm, noSeedlings, dcSeedlings, overwrite =
    FALSE)
}
\arguments{
  \item{input}{name of layer contraining number of seeds in
  soil seed bank}

  \item{disturbance}{name of disturbance layer for species}

  \item{suitability}{name of suitability layer for species}

  \item{pGermination}{proportion of seeds which will try to
  germinate, i.e. either germinate or are removed from soil
  seed bank in the process}

  \item{pGerminationDep}{\R function to determine actual
  probability of germination in dependanc of suitability
  and disturbance in cell}

  \item{seedlToDcCuts}{cuts for conversion of number of
  germinated seedlings to density classes
  \code{(parameter$res$cellSize * species$WfWdc$seedling /
  10000)}}

  \item{dcClasses}{name of density classes
  (pecies$WfWdc$dc[length(species$WfWdc$dc):1])}

  \item{pGerm}{name of number of seedlings layer}

  \item{noSeedlings}{output: name of seedlings density
  class}

  \item{dcSeedlings}{output: name of seedlings density
  class}

  \item{overwrite}{if \code{TRUE} the output file will be
  overwritten if it exists}
}
\value{
  invisible returns name of pGerm layer
}
\description{
  Basic and deterministic germination based on soil seed
  bank
}
\details{
  Calculates number of seedlings and their density class
  based on probability that seeds in seed bank germinate
  and a function
}
\author{
  Dr Rainer M Krug <rainer@krugs.de>
}

