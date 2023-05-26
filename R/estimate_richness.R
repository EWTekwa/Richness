# a wrapper function allowing mutliple data input types and choice of wheter to output point estimes or bootstrapped estimates
# created by Matt Whalen
# last update 27 Jan 2023
estimate_richness <- function( Community, boot = F, numBoot = 100, meanStates = F, Apx_detectP_terms = F ){
  if( any(class(Community) %in% c("data.frame","matrix")) ){ # note that mean states and approximated detection probability terms are only currently available when estimating single Communities (i.e., no time series or comparison of independent Communities across space)
    if( boot == T ){
      est = bootstrapRichness( Community, numBoot = numBoot )
    } else if( meanStates == T & Apx_detectP_terms == T ){
        est = RichnessEsts( Community )
      } else if( meanStates == T & Apx_detectP_terms == F ){
        est = RichnessEsts( Community )[1:2]
      } else if( meanStates == F & Apx_detectP_terms == T ){
        est = RichnessEsts( Community )[c(1:3)]
      } else if( meanStates == F & Apx_detectP_terms == F ){
        est = RichnessEsts( Community )[[1]]
      }
  } else if( class(Community) == "list" ){ # each element of the list should be a Community with taxa as columns and samples as rows
    if( boot == F ){
      estall = do.call( rbind,lapply( Community, RichnessEsts ) )
      est = do.call( rbind, estall[,1] )
    } else {
      est = do.call( rbind,lapply( Community, bootstrapRichness, numBoot = numBoot ) )
      estlong = est %>%
        mutate( id = gl( n = length(Community), k = numBoot, labels = rownames(est) ) ) %>%
        pivot_longer( Richness_raw:JK_i, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else  if( any(class(Community) == "array" & length(dim(Community)) > 2) ){ # arrays should be arrange with taxa as columns, samples as rows, and time/space in the 3rd dimension
    if( boot == F ){
      estall = do.call( rbind, apply( Community,3, RichnessEsts ) )
      est = do.call( rbind, estall[,1] )
      estlong <- est %>%
        mutate( id = as.numeric(gl( n = dim(Community)[3], k = 1 )) ) %>%
        pivot_longer( Richness_raw:JK_i, names_to = "estimate", values_to = "richness" )
      est = estlong
    } else {
      est = do.call( rbind,apply( Community, 3, bootstrapRichness, numBoot = numBoot ) )
      estlong = est %>%
        mutate( id = gl( n = dim(Community)[3], k = numBoot, labels = rownames(est) ) ) %>%
        pivot_longer( Richness_raw:JK_i, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else if( !(any(class(Community) %in% c("data.frame","matrix","list","array"))) ){
    est = c("sorry, this Community data does not appear to be in a supported format",
            "try coercing the data to a matrix, data.frame, list, or array with species as columns and spatial sampling units as rows" )
  }
      return(est)
}
