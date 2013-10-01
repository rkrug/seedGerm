## [[file:~/Documents/Projects/AlienManagementDrakensberg/sim/packages/seedGerm/seedGerm.org::*seedGermBasicGRASS][seedGermBasicGRASS:1]]
##' Basic and deterministic germination based on soil seed bank
##'
##' Calculates number of seedlings and their density class based on
##' probability that seeds in seed bank germinate and a function
##'
##' @usage seedGermBasicGRASS(input, disturbance, suitability,
##' pGermination, pGerminationDep, seedlToDcCuts, dcClasses, pGerm,
##' noSeedlings, dcSeedlings, overwrite = FALSE)
##' @name seedGermBasicGRASS
##' @title Basic seed germination
##' 
##' @param input name of layer contraining number of seeds in soil
##' seed bank
##' @param disturbance name of disturbance layer for species
##' @param suitability name of suitability layer for species
##' @param pGermination proportion of seeds which will try to
##' germinate, i.e. either germinate or are removed from soil seed
##' bank in the process
##' @param pGerminationDep \R function to determine actual probability
##' of germination in dependanc of suitability and disturbance in cell
##' 
##' @param seedlToDcCuts cuts for conversion of number of germinated
##' seedlings to density classes \code{(parameter$res$cellSize *
##' species$WfWdc$seedling / 10000)}
##' @param dcClasses name of density classes
##' (pecies$WfWdc$dc[length(species$WfWdc$dc):1])
##' @param pGerm name of number of seedlings layer
##' @param noSeedlings output: name of seedlings density class
##' @param dcSeedlings output: name of seedlings density class
##' @param overwrite if \code{TRUE} the output file will be overwritten if it exists
##' 
##' @return invisible returns name of pGerm layer
##' @author Dr Rainer M Krug <rainer@@krugs.de>
##' @export
seedGermBasicGRASS <- function(
    input,                       # layerName(species, type="seedsSoil", year=parameter$year)
    disturbance,                 # layerName(species, type="disturbances", year=parameter$year)
    suitability,                 # suitName(species)
    pGermination,                # species$pGermination
    pGerminationDep,             # species$pGerminationDep
    
    seedlToDcCuts,               # parameter$res$cellSize * species$WfWdc$seedling / 10000
    dcClasses,                   # species$WfWdc$dc[length(species$WfWdc$dc):1]
    
    pGerm,                       # layerName(species, type="pGerm", parameter$year)
    noSeedlings,                 # layerName(species, type="noSeedlings", parameter$year)
    dcSeedlings,                 # layerName(species, type="dcSeedlings", parameter$year)
    overwrite = FALSE
    ) {
    if (
        sum(
            c(pGerm, dcSeedlings, noSeedlings) %in% execGRASS("g.mlist", type="rast", mapset=".", intern=TRUE)
            ) > 0
        & !overwrite
        ) {
        stop(paste("One of the output layers exists! Please specify 'overwrite=TRUE' or use different pGerm name!"))
    } 
    germ <- readRAST6(
        c(
            sS <- input,
            dist <- disturbance,
            suit <- suitability
            ),
        NODATA=-1
        )
    
    ## calculate number of seeds which will germinate
    seedsGerm <- germ[[sS]] * pGermination
    
    ## adjust soil seed bank
    germ[[sS]] <- germ[[sS]] - seedsGerm
    
    ## successfull germination probability is pGerminationDep for the species
    germ[[4]] <- pGerminationDep(germ[[suit]], germ[[dist]])
    mode(germ[[4]]) <- "double"
    ## germ@proj4string <- parameter$proj4string
    writeRAST6(
        germ,
        vname = pGerm,
        zcol = 4,
        NODATA = -1
        )
    pGerm <- germ[[4]] 
    
    ## initialise seedl with -999
    seedl <- rep(-999, length(seedsGerm))
    
    ## set seedl to NA for all cells which have pGerm <= 0 or seeds <= 0
    seedl[pGerm <= 0] <- NA
    seedl[seedsGerm <= 0] <- NA
    
    ## set seedl to NA for all cells which have is.na(pGerm) or is.na(seedsGerm)
    seedl[is.na(seedsGerm)] <- NA
    seedl[is.na(pGerm)] <- NA
    
    ## set seedl to seedsGerm for all cells which have pGerm == 1
    rb <- (pGerm == 1)
    rb[is.na(rb)] <- FALSE
    seedl[rb] <- seedsGerm[rb]  
    
    ## identify all remaining cells (i.e. seedl == -999) cells with
    ## suitabilities 0 < suit < 1 and calculate recruits in these
    ## cells
    rb <- (seedl == -999)
    rb[is.na(rb)] <- FALSE
    oldWarn <- options()$warn
    options(warn=-1)
    seedl[rb] <- rbinom(
        n    = sum(rb),
        size = seedsGerm[rb],
        prob = pGerm[rb]
        )
    options(warn=oldWarn)
    
    seedl[seedl == 0] <- NA
    germ[[5]] <- seedl  
    mode(germ[[5]]) <- "double"
    ## germ@proj4string <- parameter$proj4string
    writeRAST6(
        germ,
        vname = noSeedlings,
        zcol = 5,
        NODATA = -1
        )
    
    ## Reclassification of noSeedlings into dcSeedlings
    dcSeedl <- as.numeric(
        as.character(
            cut(
                x = seedl,
                breaks = c(Inf, seedlToDcCuts),
                labels = dcClasses,
                include.lowest = TRUE
                )
            )
        )
    dcSeedl[dcSeedl == 0] <- NA
    germ[[5]] <- dcSeedl  
    mode(germ[[5]]) <- "double"
    ## germ@proj4string <- parameter$proj4string
    writeRAST6(
        germ,
        vname = dcSeedlings,
        zcol = 5,
        NODATA = -1
        )
}
## seedGermBasicGRASS:1 ends here

## Local Variables:
## org-babel-tangled-file: t
## buffer-read-only: t
## eval:: (auto-revert-mode)
## End:
