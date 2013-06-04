## [[file:~/Documents/Projects/AlienManagementDrakensberg/sim/packages/seedGerm/seedGerm.org::*seedGermBasic][seedGermBasic:1]]
##' Basic and deterministic germination based on seed bank
##'
##' Calculates number of seedlings and their density class based on probability that
##' seeds in seed bank germinate and a function 
##' @title 
##' @param seedsSoilLayer name of layer contraining number of seeds in soil seed bank
##' @param disturbancesLayer name of disturbance layer for species
##' @param suitabilityLayer name of suitability layer for species
##' @param pGermination proportion of seeds which will try to germinate, i.e. either germinate or are removed from soil seed bank in the process
##' @param pGerminationDep function to determine actual probability of germination in dependanc of suitability and disturbance in cell
##' @param seedlToDcCuts cuts for conversion of number of germinated seedlings to density classes (parameter$res$cellSize * species$WfWdc$seedling / 10000)
##' @param dcClasses name of density classes (pecies$WfWdc$dc[length(species$WfWdc$dc):1])
##' @param pGermLayer output: name of layer containing actual probabilities of germination per cell
##' @param noSeedlingsLayer output: name of number of seedlings layer
##' @param dcSeedlingsLayer output: name of seedlings density class (dc)
##' @return 
##' @author Dr Rainer M Krug
seedGermBasic <- function(
    seedsSoilLayer,                       # layerName(species, type="seedsSoil", year=parameter$year)
    disturbancesLayer,                    # layerName(species, type="disturbances", year=parameter$year)
    suitabilityLayer,                     # suitName(species)
    pGermination,                         # species$pGermination
    pGerminationDep,                      # species$pGerminationDep

    seedlToDcCuts,                        # parameter$res$cellSize * species$WfWdc$seedling / 10000
    dcClasses,                            # species$WfWdc$dc[length(species$WfWdc$dc):1]

    pGermLayer,                           # layerName(species, type="pGerm", parameter$year)
    noSeedlingsLayer,                     # layerName(species, type="noSeedlings", parameter$year)
    dcSeedlingsLayer                      # layerName(species, type="dcSeedlings", parameter$year) 
    ) {
    
    germ <- readRAST6(
        c(
            sS <- seedsSoilLayer,
            dist <- disturbancesLayer,
            suit <- suitabilityLayer
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
    pGerm <- germ[[4]] 
    ## germ@proj4string <- parameter$proj4string
    writeRAST6(
        germ,
        vname = pGermLayer,
        zcol = 4,
        NODATA = -1
        )
    
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
    
    ## identify all remaining cells (i.e. seedl == -999) cells with suitabilities 0 < suit < 1 and calculate recruits in these cells
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
        vname = noSeedlingsLayer,
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
        vname = dcSeedlingsLayer,
        zcol = 5,
        NODATA = -1
        )
}
## seedGermBasic:1 ends here
