#' OLD version for the main function. Does not use c++ functions. No
#' \code{fixation_at} parameter. Accepts \code{growth_rate}.
#'
#' More importantly, this model is more directly based on the UNTB and has
#' substitutions
#'
#'
#' This function simulates the changes in abundances for a given initial
#' abundance data. The model is based on the UNTB. Each iteration, individuals
#' die and are randomly substituted for new ones (without immigration or
#' speciation).

#' Additionally, the model has a dilution and growth system. After each
#' dilution, a random organism is duplicated. This happens, by default, until
#' the total abundance reaches the value previous to the dilution. But it can
#' be changed with the option abun_total
#'
#' @param counts_data
#' @param dilution
#' @param no_of_dil
#' @param abun_total Total abundance for the community to grow in each dilution-growth
#' cycle. When this value is reached, the cycle ends and another dilution (if needed)
#' happens. If NULL, \code{abun_total} is set to the initial total abundance
#' @param grow_step Number of individuals that grow each timestep (Default: 1)
#' @param keep_all_timesteps
#' @param N Substitutions happen N times each before each dilution, N being
#' by default the initial total abundance
#' @param subs_step Number of individuals substituted each time step
#'
#' @importFrom untb as.count isolate
#' @importFrom stats setNames

#' @return
#' @export
simulate_timeseries_with_sub <- function (counts_data,
                                          dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                          no_of_dil=12,
                                          N=NULL, # FIX empiezan a pasar desde justo después de la dilución (me parece lo más lógico), pero podria moverlo a justo después de alcanzar la abundancia anterior...
                                          abun_total=NULL,
                                          grow_step=1,
                                          subs_step=1,
                                          keep_all_timesteps=FALSE){
  start <- as.count(.my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)){ # establish default value.
    abun_total <- sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }
  # initialize data.frame with starting abundances
  if (keep_all_timesteps){
    trajectory <- matrix(NA, ncol=len, nrow = no_of_dil+1)
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- names(start)
    trajectory["0",]=start}

  # y el "zero counter"
  empty <- as.count(1:len)                                                       #0
  for (ind in 1:len){empty[as.character(ind)]=0}                                 #0
  names(empty) <- names(start)  #renaming is important for this function only    #0

  this_timestep <- start

  for (i in 1:no_of_dil){
    # Hago la dilución (diluimos el inóculo inicial también, como en el experimento original)
    message(paste0("before dilution: ", sum(this_timestep)))
    this_timestep <- isolate(this_timestep, (sum(this_timestep)*dilution))         # count/int ->count
    message(paste0("after dilution: ", sum(this_timestep)))
    # indicamos los ceros que hayan aparecido
    temp <- empty                                                                #0
    temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
    this_timestep <- temp[order(names(temp))]                                    #0

    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    subst_done <- 0
    if (is.null(N)) {N=abun_total}
    while ((sum(this_timestep) < abun_total) || (subst_done < N)) {
      # hago que un bicho se duplique al azar mediante sampling
      # NOTA: como estoy cogiendo del as.census, NO va a aparecer ninguno nuevo
      # de los que ya tengan abundancia 0

      ## PREPARE STEP
      if (sum(this_timestep) < abun_total) { # si aún debe seguir creciendo
        # Step definition:
        if (sum(this_timestep)==0) {
          # feb 13 2022: step to 1, and pick a random single bug anyway
          message("WARNING: less than 0.5 bugs were left after diluting! Consider changing your dilution factor.")
          message("Picking one random bug and running simulation anyway...")
          step <- 1
        } else if (sum(this_timestep) <= 3) {
          message("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.")
          message("Running simulation anyway...")
          ### Changed feb 1 2022: i had a "skip simuls for this case" but I changed
          ###                     it so it's "change step to one and throw WARNING"
          step <- 1
        } else if ((sum(this_timestep) < grow_step/2) || (sum(this_timestep) < grow_step)) { # the step is too big, reduce it
          step <- max(trunc(sum(this_timestep)/2))
          warning(paste("Step reduced to",step))
        } else {
          step <- grow_step
        }
        if ((sum(this_timestep)+step) > abun_total) { #avoid growing too much
          step=(abun_total-sum(this_timestep))}

        new_bugs <-as.count(sample(names(this_timestep),                       # count -> count
                                   size = step,
                                   prob = this_timestep,
                                   replace = TRUE))[names(this_timestep)]
        new_bugs <- new_bugs[!is.na(new_bugs)]

        new_summed=this_timestep[names(new_bugs)]+new_bugs
        this_timestep[names(new_bugs)]<-new_summed #añado al que ha nacido al censo

      }
      ## PREPARE STEP
      if (subst_done < N) { # si aun deben seguir las sustituciones.
        if ((sum(this_timestep) < subs_step/2) || (sum(this_timestep) < subs_step)) { # the step is too big, reduce it
          sstep <- max(trunc(sum(this_timestep)/2))
          # warning(paste("Step reduced to",step))
        } else {
          sstep <- subs_step
        }
        if ((sum(this_timestep)+sstep) > abun_total) { #avoid growing too much
          sstep=(abun_total-sum(this_timestep))}

        dying_bug <-as.count(sample(names(this_timestep),                       # count -> count
                                    size = sstep,
                                    prob = this_timestep,
                                    replace = TRUE))[names(this_timestep)]
        dying_bug <- dying_bug[!is.na(dying_bug)]

        without_dead=this_timestep[names(dying_bug)]-dying_bug
        this_timestep[names(dying_bug)]<-without_dead # kill it

        temp <- empty; temp[names(this_timestep)] <- this_timestep #fix zeros    #0

        new_bugs <-as.count(sample(names(this_timestep),                       # count -> count
                                   size = sstep,
                                   prob = this_timestep,
                                   replace = TRUE))[names(this_timestep)]
        new_bugs <- new_bugs[!is.na(new_bugs)]
        new_summed=this_timestep[names(new_bugs)]+new_bugs
        this_timestep[names(new_bugs)]<-new_summed #añado al que ha nacido al censo

        subst_done <- subst_done + 1
      }
    }

    # una vez crecidos, se puede diluir de nuevo: se repite el bucle

    # pero primero "secuenciamos" (guardamos las abundancias actualizadas)
    if (keep_all_timesteps){
      trajectory[as.character(i),]=this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep%>%setNames(names(start)))
  }
}
