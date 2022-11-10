sourceCpp("./growth.cpp")

create_counts <- function (exp, # abundances
                           map, # metadata about m_inic, replicates, sample names...
                           m_inic, # original sample names
                           allow_empty_transfers=FALSE,
                           select_transfers=NULL, # vector of transfers to select
                           orig="ORIG", # original sample tag column name in map
                           transfer="T",# transfer tag column name in map
                           sa="SA"     # sample name tag column name in map
){ 
  # This function parses and puts together in a list of data.frames all the 
  # count data for time series abundance data. For each timestep / transfer, it
  # saves the mean value of all available replicates.
  
  # If there are no available replicates for a given transfer, it saves the data
  # for the previous transfer (in other words, assumes no abundance changes took
  # place). This assumption can be overridden with the option 
  # "allow_empty_transfers", that will set "NA" instead.
  
  # If select_transfers is not NULL, is has to be a vector. If you want to pick
  # only transfers 0, 2 and 4, you must input c(0,2,4)
  
  counts_ <- list()
  if (!is.null(select_transfers)){
    transfers=select_transfers
  } else {
    transfers=0:length(unique(map[[transfer]]))      # includes transfer "0" (original)
  }
  for (i in m_inic) {
    counts_[[i]]<-data.frame(matrix(ncol = 0, nrow = nrow(exp))) # prepare empty DF
    for (t in transfers){
      ### separo por muestra inicial los datos finales también
      temp <- exp[map[map[orig]==i & map[transfer]==t,][,sa]]
      if (ncol(temp)==0){  # if we don't have data for this m_inic in this transfer, 
        # skip and repeat the data for the last transfer instead.
        if (!allow_empty_transfers){
          counts_[[i]]=cbind(counts_[[i]],
                             setNames(counts_[[i]][ncol(counts_[[i]])],nm=t)) 
        } else {
          counts_[[i]]=cbind(counts_[[i]],
                             setNames(as.data.frame(rep(NA,nrow(counts_[[i]]))),nm=t))
        }
      } else {
        ### calculate the mean for all replicates for that transfer
        r <- rowSums(temp)/dim(temp)[2] # might return float
        ### finally, I add the mean data for that transfer to the final counts_ list
        counts_[[i]]=cbind(counts_[[i]],
                           setNames(as.data.frame(r),nm=t)) #añado la nueva transfer
        counts_[[i]][is.na(counts_[[i]])]=0
      }
    }                                                                                   
  }  
  return(counts_)# FIX en el script usar create counts para parsear el resultado de la simulación, guardado ya en dataframe...
}




simulate_timeseries <- function (counts_data,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 fixation_at=1,
                                 abun_total=NULL,
                                 grow_step=1,
                                 keep_all_timesteps=FALSE,
                                 force_continue=FALSE) { 
  
  start <- as.count(my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent
  
  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)){ # establish default value.
    abun_total <- sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }
  
  # if keep_all_timesteps, initialize data.frame with starting abundances
  if (keep_all_timesteps){
    trajectory <- matrix(NA, ncol=len, nrow = no_of_dil+1)
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- names(start)
    trajectory["0",]=start
  } 
  
  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0
  
  this_timestep <- start
  dil <- 0
  while (dil < no_of_dil &
         (max(this_timestep)/sum(this_timestep) < fixation_at)){
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep)*dilution)==0) {
      if (force_continue) {
        message("WARNING: less than 1 bug left after diluting! Consider changing your dilution factor.")
        message("Picking one random bug and running simulation anyway...")   
        i <- c(1:length(this_timestep))[.Internal(sample(
          x = length(this_timestep),
          size = 1,
          replace = FALSE,
          prob = this_timestep))]
        this_timestep    <- empty
        this_timestep[i] <- this_timestep[i] + 1
      } else {
        write("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.", stdout())
        stop("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.")
      }
      ## If not CASE 1, dilute normally
    } else {
      ## CASE 2 -- throw a warning if not many bugs left
      if (sum(this_timestep) <= 3) {
        write("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.", stdout())  # TODO 3 es ridículo, pensar uno más grande
        write("Running simulation anyway...", stdout())
      }
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep)*dilution,
        replace = TRUE,
        prob = this_timestep))]) 
      
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                    #0
    }
    
    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta 
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    ns <- names(this_timestep)
    this_timestep <- as.vector(this_timestep)
    while (sum(this_timestep) < abun_total) {
      this_timestep <- growth(this_timestep, abun_total, grow_step)                # Rcpp function
    }
    names(this_timestep) <- ns
    dil <- dil + 1
    
    if (keep_all_timesteps){
    # una vez crecidos, se puede diluir de nuevo: siguiente iteración del bucle
    # pero si keep_all(...), antes "secuenciamos" (guardamos las abundancias)
      trajectory[as.character(dil),]=this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  if (max(this_timestep)/sum(this_timestep) >= fixation_at) {
    write(paste0(fixation_at*100, "% fixation of ",
                 names(this_timestep)[this_timestep!=0], " after ", dil, " dilutions."), stdout())
  }
  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep)
  }}



simulate_timeseries_old <- function (counts_data,
                                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                     no_of_dil=12,
                                     abun_total=NULL,
                                     grow_step=1,
                                     keep_all_timesteps=FALSE,
                                     growth_rate=NULL){
  # This function simulates the changes in abundances for a given initial 
  # abundance data. The model has a dilution and growth system. After each 
  # dilution, a random organism is duplicated. This happens, by default, until 
  # the total abundance reaches the value previous to the dilution. But it can 
  # be changed with the option abun_total. "grow_step" is the number of 
  # individuals that grow each timestep.
  
  # keep_all_timesteps==FALSE saves a bit of RAM if there are many dilutions.
  # Only advised to make it TRUE for plotting.
  
  start <- as.count(my_transpose(counts_data))
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
    trajectory["0",]=start
  } 
  
  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0
  
  this_timestep <- start
  
  message(paste0("Original abundance: ", sum(this_timestep)))
  message(paste0("After dilution: ", dilution*sum(this_timestep)))
  
  for (i in 1:no_of_dil){
    # Hago la dilución (diluimos el inóculo inicial también, como en el experimento original)
    
    # new_bugs <-c(1:length(this_timestep))[.Internal(sample(
    #   x = length(this_timestep),
    #   size = sum(this_timestep)*dilution, 
    #   replace = TRUE,
    #   prob = this_timestep))]
    # 
    # for (i in new_bugs) {
    #   this_timestep[new_bugs] <- this_timestep[new_bugs] + 1 # añado al que ha nacido al censo
    # }
    
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep)*dilution)==0) {
      # feb 13 2022: step to 1, and pick a random single bug anyway
      message("WARNING: less than 1 bug left after diluting! Consider changing your dilution factor.")
      message("Picking one random bug and running simulation anyway...")   
      i <- c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = 1,
        replace = FALSE,
        prob = this_timestep))]
      this_timestep    <- empty
      this_timestep[i] <- this_timestep[i] + 1
      print(this_timestep)
      ## If not CASE 1, dilute normally
    } else {
      ## CASE 2 -- throw a warning if not many bugs left
      if (sum(this_timestep) <= 3) {
        message("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.")
        message("Running simulation anyway...")
        ### Changed feb 1 2022: i had a "skip simuls for this case" but I changed
        ###                     it so it's "change step to one and throw WARNING"
      }
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep)*dilution,
        replace = TRUE,
        prob = this_timestep))])
      # indicamos los ceros que hayan aparecido
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                    #0
    }
    
    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta 
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    while (sum(this_timestep) < abun_total) {# si aún debe seguir creciendo
      if (sum(this_timestep) < grow_step) { # the step is too big, reduce it
        step <- max(trunc(sum(this_timestep)/2))
        # warning(paste("Step reduced to",step))
      } else {
        step <- grow_step
      }
      
      ## Avoid growing too much (when step>1)
      if ((sum(this_timestep)+step) > abun_total) {
        step=(abun_total-sum(this_timestep))}
      
      ## Once the step is decided, sample
      new_bugs <- c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = step, 
        replace = TRUE,
        prob = this_timestep))]
      for (i in new_bugs) {
        this_timestep[i] <- this_timestep[i] + 1
      } # añado al que ha nacido al censo
      
    }
    
    # una vez crecidos, se puede diluir de nuevo: se repite el bucle
    # pero primero "secuenciamos" (guardamos las abundancias actualizadas)
    if (keep_all_timesteps){
      trajectory[as.character(i),]=this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  message(paste0("Reached abundance of: ", sum(this_timestep)))
  
  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep%>%setNames(names(start)))
  }
}



# FIX: translate comments
simulate_timeseries_with_sub <- function (counts_data,
                                          dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                          no_of_dil=12,
                                          N=NULL,
                                          abun_total=NULL,
                                          grow_step=1,
                                          subs_step=1,
                                          keep_all_timesteps=FALSE){ 
  # This function simulates the changes in abundances for a given initial 
  # abundance data. The model is based on the UNTB. Each iteration, individuals 
  # die and are randomly substituted for new ones (without immigration or 
  # speciation). 
  
  # "with_sub" == Substitutions happen N times each before each dilution, N being
  # by default the initial total abundance. 
  # FIX empiezan a pasar desde justo después de la dilución (me parece lo más lógico), pero podria moverlo a justo después de alcanzar la abundancia anterior...
  
  # Additionally, the model has a dilution and growth system. After each 
  # dilution, a random organism is duplicated. This happens, by default, until 
  # the total abundance reaches the value previous to the dilution. But it can 
  # be changed with the option abun_total
  
  # "grow_step" is the number of individuals that grow each time step. subs_step
  # is the number of individuals substituted each time step
  
  # keep_all_timesteps==FALSE saves RAM. Only advised to make it TRUE for 
  # plotting. 
  start <- as.count(my_transpose(counts_data))
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






simulate_timeseries_gr <- function (counts_data,
                                    dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                    no_of_dil=12,
                                    abun_total=NULL,
                                    grow_step=1,
                                    keep_all_timesteps=FALSE,
                                    growth_rate=NULL){
  # This function simulates the changes in abundances for a given initial 
  # abundance data. The model has a dilution and growth system. After each 
  # dilution, a random organism is duplicated. This happens, by default, until 
  # the total abundance reaches the value previous to the dilution. But it can 
  # be changed with the option abun_total. "grow_step" is the number of 
  # individuals that grow each timestep.
  
  # keep_all_timesteps==FALSE saves a bit of RAM if there are many dilutions.
  # Only advised to make it TRUE for plotting.
  
  ## growth_rate :: (added oct 18 2021)
  ##   Indicates the name of an optional column that determines how much a bug 
  ##   grows. Each time it grows, it won't grow by 1 but by this factor. It
  ##   indicates the grow rate in percent mass.
  ##   WARNING: It can't be a decimal number because of how the sampling works.
  ##   That means the count_data data.table, which includes the growth rates, 
  ##   must be multiplied by a number that turns all rates into an integer. The
  ##   abun_total AND grow_step arguments must be multiplied by this number, too   
  ###     > 詳しく: supongamos growth rates (gr) de 2 decimales, como 0.23. 
  ###       Entonces nuestro number será 100. Entonces cuando el bicho con gr de
  ###       0.23 crezca 100, antes de sumarlo al total le haremos 100*0.23 = 23.
  
  
  # Take growth rates if defined
  if (!is.null(growth_rate)) {
    growths <- counts_data[growth_rate]
    counts_data <- counts_data[!names(counts_data)==growth_rate]
  }
  
  start <- as.count(my_transpose(counts_data))
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
    trajectory["0",]=start
  } 
  
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
    total_bugs_timestep <- sum(this_timestep)
    while (total_bugs_timestep < abun_total) {# si aún debe seguir creciendo
      # hago que un bicho se duplique al azar mediante sampling
      # NOTA: como estoy cogiendo del as.census, NO va a aparecer ninguno nuevo 
      # de los que ya tengan abundancia 0
      
      ## FIX HERE
      if (sum(this_timestep)==0) {
        # feb 13 2022: step to 1, and pick a random single bug anyway
        message("WARNING: less than 0.5 bugs were left after diluting! Consider changing your dilution factor.")
        message("Picking one random bug and running simulation anyway...")   
        step <- 1
        this_timestep[sample(1:length(this_timestep), 1)] <- 1
      } else if (sum(this_timestep) <= 3) {
        message("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.")
        message("Running simulation anyway...")
        ### Changed feb 1 2022: i had a "skip simuls for this case" but I changed
        ###                     it so it's "change step to one and throw WARNING"
        step <- 1
      } else if (total_bugs_timestep < grow_step) { # the step is too big, reduce it
        step <- max(trunc(total_bugs_timestep/2))
        # warning(paste("Step reduced to",step))
      } else {
        step <- grow_step
      }
      
      if ((total_bugs_timestep+step) > abun_total) { #avoid growing too much
        step=(abun_total-total_bugs_timestep)}
      
      if (!is.null(growth_rate)) {
        prob = growths[names(this_timestep),]*this_timestep
      } else {
        prob = this_timestep
      }
      
      ## Once the step is decided, sample
      new_bugs <-c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = step, 
        replace = TRUE,
        prob = this_timestep))]
      
      for (i in new_bugs) {
        this_timestep[new_bugs] <- this_timestep[new_bugs] + 1 # añado al que ha nacido al censo
      }
      
      total_bugs_timestep <- sum(this_timestep)
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




make_winner_report <- function(all_final_abund, # data.frame with the final 
                               # abundances (rows) for every
                               # simulation
                               initial_RA, # initial relative abundances
                               ignore_zeroes=TRUE) { # when calculating the mode
  # e.g.:
  # > abund_report["4456889"]
  # 4456889
  # V1        NA
  # V2        NA
  # V3         2
  # V4        NA
  # mode       2
  
  # For each row, obtain a list of the position of each otu when ordered by
  # abundance. Then, add a row with the most common position for each OTU,
  # another with the percentage or simulations where it has that position and
  # two final rows with the initial absolute abundances and the survival rate.
  #
  #   ITER     | fulano | mengano   | merengano | etc
  # -----------|--------|-----------|----------
  #   1        | 1      | 2         | 3
  # -----------|--------|-----------|----------
  #   2        | 2      | 1         | 3
  # -----------|--------|-----------|----------
  #   3        | 1      | 2         | 3
  # -----------|--------|-----------|----------
  #   4        | 1      | 2         | 3
  # -----------|--------|-----------|----------
  #   5        | 1      | 2         | 3
  # -----------|--------|-----------|----------
  # mode (most |        |           |
  # frequent   |   1    | 2         | 3        <- table is ordered by this field
  # position)  |        |           |
  # -----------|--------|-----------|----------
  # percentage |        |           |
  # of times   |   80%  | 80%       | 100 % 
  # with that  |        |           |
  # position   |        |           |
  # -----------|--------|-----------|----------
  # initial RA | mucha  | no tanta  | menos
  # -----------|--------|-----------|----------
  # survival   | mucha  | no tanta  | menos
  #   rate     |        |           | 
  
  # data.frame with positions for each otu in each simulation
  tryCatch(
    expr = {abund_report <- -all_final_abund %>%
      apply(MARGIN = 1, rank,
            ties.method="max",# Extinct species will have the biggest index possible
            simplify = FALSE) %>%# 08feb2022 avoid mistakes when there's only one column
      as.data.frame() %>%        # the rank will always be "NA" but it will be CORRECT
      my_transpose()
    }, error=function(cond){
      message("We'll just run the function without 'simplify'...") # TODO 13feb this failed when evaluating "all", maybe it fails when redundant? or the "class" is different when it's "all"?
      abund_report <- -all_final_abund %>%
        apply(MARGIN = 1, rank,
              ties.method="max") %>% # Extinct species will have the biggest index possible
        # 08feb2022 avoid mistakes when there's only one column
        as.data.frame() %>%        # the rank will always be "NA" but it will be CORRECT
        my_transpose()
    }
  )
  
  # 08feb2022 fix the row indexes !
  rownames(abund_report) <- 1:nrow(abund_report)
  
  # If we are not going to take into account zeros, we remove the biggest
  # indexes for calculating the mode (last position is almost always 0)
  if (ignore_zeroes){
    abund_report[abund_report==ncol(abund_report)] <- NA
  }
  
  # add a row with the mode of each col (custom function "mode")
  abund_report <- abund_report %>% apply(MARGIN = 2, mode) %>%
    as.data.frame() %>% 
    my_transpose() %>%
    rbind(abund_report,"mode"=.) # add mode
  
  # add percentage of times with that position
  abund_report["mode prevalence",] <- abund_report %>% 
    apply(MARGIN=2,
          FUN=function(col){
            len=length(col)
            if(!is.na(col["mode"])){ # if the mode is not NA, calculate %
              return(sum(col[1:(len-1)]==col["mode"],na.rm = T) / sum(!is.na(col[1:len-1]))* 100)
            } else { return(NA)}
          }
    )
  
  # add RA
  abund_report["RA",] <- colnames(abund_report) %>%
    initial_RA[.,] 
  
  abund_report["survival rate",] <- abund_report %>% 
    apply(.,MARGIN=2,
          FUN=function(col){         
            len=length(col)
            return((sum(!is.na(col[1:(len-3)])) / (len-3)) * 100)
          }
    )
  
  # order by mode
  abund_report <- abund_report %>% .[order(.["mode prevalence",],decreasing = T)] 
  return(abund_report)
}



PA_PD_df_maker<- function(PCG_name=NULL) {
  df <- data.frame(matrix(ncol = 4, nrow = 10),# prepare empty DF
                   row.names = c("PA_abun","PD_abun",
                                 "PA_shannon","PD_shannon",
                                 "PA_pielou","PD_pielou",
                                 "PA_obs_rich","PD_obs_rich",
                                 "PA_faith","PD_faith")) %>%
    setNames(c("all_PCGs", PCG_name))
}




PA_PD_abundances <- function(simul_abund,
                             real_abund,
                             outfile) { # names for data.frame if 
  # PA_PD_results_dataframe==NULL
  # Given a data.frame of abundata (real and simulated data), creates a 
  # metadata data.frame (groups are 0=simulated and 1=real data) and fills in 
  # a given data.frame with the PERMANOVA and PERMADISP results from comparing
  # these two groups.
  
  # If a data.frame to fill in is not provided, it well be created with 
  # PA_PD_df_maker().
  # -------------------------------------------------------
  # select common OTUs; it can happen that initial data (used for simul_abund)
  # doesn't include an OTU that appears later in the final data of real_abund.
  # This happens when the OTU was under the detectiom limit for sequencing in 
  # the initial samples.          ## glc -> init 192 to fin 167 otus.
  real_abund <- real_abund[colnames(simul_abund)]            
  # Remove NA
  simul_abund = simul_abund %>% na.omit
  real_abund  = real_abund  %>% na.omit
  
  # Check abundance of real data. If abundance is 0 for a sample, it gets 
  # deleted and the user is warned. If all samples have 0 abundance, stop.
  sums <-  rowSums(real_abund)
  sums2 <-  rowSums(simul_abund)
  if (any(sums==0)){
    message(paste0("There were ",sum(sums==0)," samples with 0 sums in the real data"))
    real_abund <- real_abund[names(sums[sums!=0]),,drop=F] # remove sum=0 samples   
    if (all(sums==0)){
      stop("real_data is empty, PERMANOVA cannot be carried out")                # FIX ?
    }}
  sums2 <-  rowSums(simul_abund)
  if (any(sums2==0)){
    message(paste0("There were ",sum(sums2==0)," samples with 0 sums in the simulated data"))
    simul_abund <- simul_abund[names(sums2[sums2!=0]),,drop=F] # remove sum=0 samples
    if (all(sums2==0)){
      stop("simul_data is empty, PERMANOVA cannot be carried out")                # FIX ? If this happens, its a good idea to just remove the "limit" earlier so we can use more simulations
    }}
  
  simul_names <- sprintf("simul%s",seq(1:nrow(simul_abund))) # rename
  rownames(simul_abund) <- simul_names
  
  abundata <- simul_abund %>% rbind(.,                        ## real data
                                    real_abund[colnames(.)])  ## simulated data
  
  metadata <- abundata %>%
    rownames() %>%
    str_detect(pattern="simul",negate=T) %>%           # FALSE if simulated, 
    as.numeric() %>%                                       # TRUE if real data 
    as.data.frame %>%
    setNames(.,"group")
  
  ## PA
  resultPA <- adonis(abundata ~ group, data=metadata, permutations = 9999)
  PA <- resultPA$aov.tab$`Pr(>F)`[1]
  
  ## PD
  abundata_dist <- abundata %>%# my_transpose() %>%
    vegdist(.,method="bray") 
  permdisp <- betadisper(abundata_dist, group = metadata$group)
  resultPD <- permutest(permdisp, permutations = 9999)
  PD <- resultPD$tab$`Pr(>F)`[1]
  
  ## plot PD
  png(paste0(output,"/",outfile,".png"))
  par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
  plot(permdisp, main="PCoA")
  boxplot(permdisp, main="Distance to centroids")  
  mtext(paste0(outfile), 
        outer=T, cex=1.2)        
  dev.off()
  
  return(list(PA,PD))
}

plot_ecdf_KS <- function(data_simul,data_real, metric, m, PCG="all") {
  simul_data = data_simul[[metric]] %>% setNames(.,nm = ".") %>% .[["."]]
  real_data  = data_real [[metric]] %>% setNames(.,nm = ".") %>% .[["."]]
  
  # fuente https://stackoverflow.com/questions/39162178/kolmogorov-smirnov-plot-in-r-ggplot
  group <- c(rep("simul_data", length(simul_data)), rep("real_data", length(real_data)))
  dat <- data.frame(KSD = c(simul_data,real_data), group = group)
  cdf1 <- ecdf(simul_data) 
  cdf2 <- ecdf(real_data) 
  
  minMax <- seq(min(simul_data, real_data, na.rm = T), max(simul_data, real_data, na.rm = T), length.out=length(simul_data)) 
  x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
  y0 <- cdf1(x0) 
  y1 <- cdf2(x0) 
  
  pva = ks.test(simul_data,real_data)$p.value %>% format(round(., 2), nsmall = 2)
  
  print(
    ggplot(dat, aes(x = KSD, group = group, colour = group, linetype=group))+
      # scale_linetype_manual(values=c("blank", "dashed"), labels=group) +
      stat_ecdf(size=1,geom="point") +
      scale_color_manual(values=c("blue","gold"), labels= c("real_data","simul_data"))+ 
      xlab(metric) +
      ylab("Cumulative Distibution") +
      annotate(geom = "text", label=pva, size=3, color="darkred",
               x = x0[1]+0.5 , y= y0[1]) +
      geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
                   linetype = "dashed", color = "coral",size=0.1) +
      geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=0.1) +
      geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=0.1) +
      ggtitle(paste0("K-S Test: Simul_data / Real_data -- ",metric," ",m," PCG: ",PCG))
  )
  # for (d in c(1,2)){
  #   data=list(simul_data,real_data)[[d]]
  #   name=c("simul","real")[d]
  #   plots[[name]] <-
  #     `stat_ecdf`(mapping = aes_string(data[["."]],
  #                                                   color=as.factor(name)),
  #                              geom = "point")
  # }
  # 
  # print(ggplot()+ #empty backbone
  #         plots + #datos
  #         labs(title=paste0("Empirical Cumulative Density Function -- ",metric),
  #              y = paste0("(%)"), x=metric) +
  #         scale_fill_manual("",
  #                            breaks=names,
  #                            values=colors) +
  #         theme(legend.position="bottom")
  # )
}

plot_ecdf_KS2 <- function(data_simul,data_real, metric, m, PCG="all") {
  simul_data = data_simul[[metric]] %>% setNames(.,nm = ".") %>% .[["."]]
  real_data  = data_real [[metric]] %>% setNames(.,nm = ".") %>% .[["."]]
  
  # fuente https://stackoverflow.com/questions/39162178/kolmogorov-smirnov-plot-in-r-ggplot
  group <- c(rep("simul_data", length(simul_data)), rep("real_data", length(real_data)))
  dat <- data.frame(KSD = c(simul_data,real_data), group = group)
  cdf1 <- ecdf(simul_data) 
  cdf2 <- ecdf(real_data) 
  
  minMax <- seq(min(simul_data, real_data, na.rm = T), max(simul_data, real_data, na.rm = T), length.out=length(simul_data)) 
  x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
  y0 <- cdf1(x0) 
  y1 <- cdf2(x0) 
  
  pva = ks.test(simul_data,real_data)$p.value %>% format(round(., 2), nsmall = 2)
  
  print(
    ggplot(dat, aes(x = KSD, group = group, colour = group, linetype=group))+
      # scale_linetype_manual(values=c("blank", "dashed"), labels=group) +
      # ECDF: ----
    stat_ecdf(size=1,geom="point") +
      scale_color_manual(values=c("blue","gold"), labels= c("real_data","simul_data"))+ 
      xlab(metric) +
      ylab("Cumulative Distribution") +
      # ----------
    
    # box: -----
    geom_boxplot(outlier.colour="red", 
                 color="black",
                 alpha=0.5,
                 fill=c("blue","gold"),
                 outlier.shape=16,
                 outlier.size=2, notch=FALSE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) +
      # ----------
    
    # KS: ------
    annotate(geom = "text", label=pva, size=3, color="darkred",
             x = x0[1]+0.5 , y= y0[1]) +
      geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
                   linetype = "dashed", color = "coral",size=0.1) +
      geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=0.1) +
      geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=0.1) +
      # ----------
    
    # title: ---
    ggtitle(paste0("K-S Test: Simul_data / Real_data -- ",metric," ",m," PCG: ",PCG))
    # ----------
  )
  # for (d in c(1,2)){
  #   data=list(simul_data,real_data)[[d]]
  #   name=c("simul","real")[d]
  #   plots[[name]] <-
  #     `stat_ecdf`(mapping = aes_string(data[["."]],
  #                                                   color=as.factor(name)),
  #                              geom = "point")
  # }
  # 
  # print(ggplot()+ #empty backbone
  #         plots + #datos
  #         labs(title=paste0("Empirical Cumulative Density Function -- ",metric),
  #              y = paste0("(%)"), x=metric) +
  #         scale_fill_manual("",
  #                            breaks=names,
  #                            values=colors) +
  #         theme(legend.position="bottom")
  # )
}


plot_boxplot <- function(data_simul,data_real, metric) {
  bpdata = rbind(data_real[[metric]] %>% cbind(.,"tag"=rep("real",ncol(.))),
                 data_simul[[metric]] %>% cbind(.,"tag"=rep("simul",ncol(.)))) %>% setNames(c("val","tag"))
  
  print(
    ggplot(bpdata, aes(x = tag, y=val, fill=tag))+
      geom_boxplot(outlier.colour="red", color="black",
                   outlier.shape=16,
                   outlier.size=2, notch=FALSE) +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)
  )
}

# 24 sep 2021 gamma_dist_model.R
create_gamma_distr_simuls <- function(data, n) {
  ## data --> real (or not) abundance datasets from where we will obtain
  ##          a mean abundance and its variance for every present OTU
  ## n -----> number of simulations wanted
  data["mean"] = apply(data,MARGIN = 1,FUN = mean) %>% as.data.frame #x_
  data["var"]  = apply(data,MARGIN = 1,FUN = var)  %>% as.data.frame #o2
  beta = data["mean"]**2/data["var"]
  beta[is.na(beta)]<-0
  data["shape"]=beta #a
  data["scale"]=1/(beta/data["mean"]) #s
  data["scale"][is.na(data["scale"])] <-0
  
  
  simul <-apply(data[,c("shape","scale")], MARGIN = 1, FUN=function(row){
    rgamma(n=n, shape=row[1] , scale=row[2])})
}
# 30-sep-2021 2021-09-22_gamma_distr_model n2_correlations.R
# método de Cordero de usar simulaciones con distribuciones gamma para
# dar significancia a correlaciones entre OTUs y hacer buenos grafos
# ++++++++++++++++++++++++++++
# flattenCorrMatrix http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat=NULL) { # incluyo opcion de quitar pmat
  ut <- upper.tri(cormat)
  if (is.null(pmat)) {
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =cormat[ut]
    )
  } else {
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =cormat[ut],
      p = pmat[ut]
    )
  }
}


# aplicable a una pmat tmb
flattenCorrMatrixCor <- function(cormat_list) {
  first <- flattenCorrMatrix(cormat_list[[1]])
  next_ones <- lapply(cormat_list[-1], FUN = function(cormat){
    ut <- upper.tri(cormat)
    data.frame(cor=cormat[ut])
  })
  cbind(first,next_ones)%>%setNames(c("row","column",c(names(cormat_list))))
}


unflatten <- function(orig,colnam=c('otu1','otu2'),values="weight") {
  cn1 = colnam[1]
  cn2 = colnam[2]
  
  df=data.frame(row.names = c(orig[[cn1]],orig[[cn2]])%>%unique)
  for (o1 in unique(orig[[cn1]])) {
    for (o2 in unique(orig[[cn2]])) {
      tryCatch(expr={
        df[o1,o2] = orig[orig[[cn1]]==o1 & orig[[cn2]]==o2,values]
        df[o2,o1] = orig[orig[[cn1]]==o1 & orig[[cn2]]==o2,values]
      },error=function(cond){
        # message(cond)
        # message(o1)
        df[o1,o2] = 0
        df[o2,o1] = 0
      }
      )
    }
  }
  df[is.na(df)] <- 0
  df <-df[colnames(df),colnames(df)]
  return(df)
}
