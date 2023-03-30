start_time <- Sys.time()
## =============================================================================
##                                SIMULACIÓN
## =============================================================================
## El UNTB consiste en una serie de comunidades locales neutrales que reciben
## inmigración de una metacomunidad también neutral. El modelo original
## sustituye, en cada generación/timestep, DE MEDIA una vez a cada organismo.
## Esto significa que en cada timestep habrá N iteraciones, siendo N el número
## total de bichos, en las que se escogerá al azar un individuo que morirá y
## será sustituido por un nuevo individuo de una especie ya presente (división
## celular) o bien de una no presente en la comunidad local pero sí en la
## metacomunidad (inmigración).

## En mi modelo hay un par de cambios: en primer lugar se trata de comunidades
## aisladas, no hay metacomunidad ni inmigración (son aisladas porque estoy
## replicando experimentos de laboratorio); en segundo lugar, cada vez que se
## alcanza cierta abundancia total hago una dilución, simulando los transfers
## del experimento del wet-lab. Esto quiere decir que no hago N iteraciones o
## sustituciones, sino tantas como haga falta para alcanzar la misma abundancia
## total que había antes de la dilución. Hay 12 transfers, contando el inóculo
## inicial.

## Las abundancias relativas en cada timestep se guardan en una matriz y
## posteriormente se hace la media para cada OTU en cada momento y se calculan
## la richness, evenness y diversidad, medidas de las cuales también se calcula
## una media.

## =============================================================================
## "disclaimer"
## =============================================================================
## Hacemos 2 asunciones que hay que entender bien:
##
##   1) Asumimos que tras cada transfer se ha alcanzado ya la proporción final
##   conocida de cada PCG. Obviamente en la vida real no tiene por qué ser así,
##   por ejemplo: que tras 7 transfers haya 75% de Enterobacteriaceae y 25% de
##   Pseudomonadaceae pero que tras la 3ª haya todavía 55% y 45%.
##   [En Talavera et al. 2022 se tuvo en cuenta simulando a partir de una
##   comunidad inicial ya estable]
##
##   2) Asumimos que es lo mismo diluir un cultivo de 10000 que son un surtido
##   de varias PCGs que diluir esas PCGs por separado, cuando en la vida real no
##   es exactamente así: en la vida real habrá ruido. Por ejemplo, si tengo mi
##   cultivo de 100000 y 1000 son Acetobacter y quiero diluir un 20%, en esta
##   pipeline siempre serán seleccionados 200 de Acetobacter. En la vida real,
##   sin embargo, a veces serán 200, a veces 190, a veces 211... Ese ruido se
##   pierde. Dicho esto, teniendo esto en cuenta se detiene el script: si el
##   número de individuos de un PCG es demasiado bajo, se avisa al usuario de
##   que el factor de dilución debe cambiarse para evitar una posible extinción
##   del mismo, que sí podría ocurrir en la vida real.
## =============================================================================

# ===========
# mis datos [glc]
# ===========
library("tidyverse")
library("dilgrowth")
library("gsubfn")
library("parallel")
library("optparse")
library("dplyr") # para bind_rows tras paralelización
option_list <- list(
  # input params
  make_option(c("-a", "--abuntable"), type = "character",
              default = NULL,
              help = "Abundance table"),
  make_option(c("--skip_lines_abuntable"), type = "integer",
              default = 0,
              help = "How many lines of abuntable should we skip (read.csv()'s 'skip'). Should be 1 for Qiime 1 output, for instance"),
  make_option(c("-s", "--sample"), type = "character",
              default = NULL,
              help = "Sample name in the abundance table"),
  make_option(c("--subset"), type = "character",
              default = NULL,
              help = "Colon-separated list of OTUs to be selected from the abuntable (e. g. those belonging to a given PCG). NULL by default, meaning no selection)"),
  # simul params
  make_option(c("--dilution"), type = "character",
              default = "8 * 10 ** (-3)",
              help = "Dilution before each simulated transfer. Can be an equation (e.g. 10**(-3))"),
  make_option(c("--no_of_dil"), type = "integer",
              default = 12,
              help = "Number of dilutions/transfers"),
  make_option(c("--fixation_at"), type = "double",
              default = 1,
              help = "The simulations stop when only one bug is fixed at 100% [fixation_at=1], since it's meaningless to keep going when only one species is left. Nevertheless, we can consider an OTU is fixed when it has reached a lower relative abundance (e.g. fixation_at=0.95) and stop earlier."),
  make_option(c("--save_all"), type = "logical",
              default = FALSE,
              help = "save all intermediate states of the simulations? default FALSE"),
  make_option(c("--grow_step"), type = "double",
              default = 1,
              help = "How many bugs are born on each iteration. 1 by default"),
  make_option(c("--is_grow_step_a_perc"), type = "logical",
              default = FALSE,
              help = "If FALSE, grow_step is taken as a fixed value, so the step will always be the same. If TRUE, it's read as a percentage - the step will be changed proportionally to the community size. If grow_step is 0.02, 2% of the members in the community will grow the next iteration. FALSE by default"),
  make_option(c("--fix_percentage"), type = "logical",
              default = TRUE,
              help = "If TRUE, uses a fixed percentage (relative abundance, given by --perc). If FALSE, randomly chooses a percentage in each simulation. These percentages are taken from a map given by --perc_map. Default: TRUE"),
  make_option(c("--perc"), type = "character",
              default = "1",
              help = "Desired final PCG abundance (if --fix_percentages==TRUE). 1 by default."),
  make_option(c("--perc_map"), type = "character",
              default = NULL,
              help = "Percentage map, see --fix_percentage"),
  make_option(c("--cores"), type = "integer",
              default = 1,
              help = "Number of cores to use in parallelization processes (mclapply). Default: 4.",
              metavar = "integer"),
  # output params
  make_option(c("--no_of_simulations"), type = "integer",
              default = 1,
              help = "Number of simulations"),
  make_option(c("--outputname"), type = "character",
              default = "out",
              help = "name for output .csv file"),
  make_option(c("--outdir"), type = "character",
              default = "results",
              help = "output directory")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

abuntable  <- opt$abuntable  # e.g. "./original_100percArbol/Tree/0.99/table.from_biom_0.99.txt"
skip_lines_abuntable <- opt$skip_lines_abuntable
s <- opt$sample      # e.g. "X2", col name from map
subset_otus <- opt$subset

outdir <- opt$outdir # e.g. "my_neutral_model_v2_test_16_simuls_8_cores"
outputname <- opt$outputname # e.g. "X2_rep4"
cores <- opt$cores  # e.g. 16

save_all  <- opt$save_all
grow_step <- opt$grow_step
is_grow_step_a_perc <- opt$is_grow_step_a_perc
if (is_grow_step_a_perc && (grow_step == 1)) {
  message("WARNING: grow_step has been defined as a percentage of value 1. All organisms will duplicate every iteration. If you meant for grow_step to be a fixed value, set is_grow_step_a_perc to FALSE instead.")
}

dilution <-  eval(parse(text=opt$dilution)) # e.g. 8 * 10 ** (-3)
no_of_dil <- opt$no_of_dil
fixation_at <- opt$fixation_at

no_of_simulations <- opt$no_of_simulations

fix_percentage <- opt$fix_percentage
if (fix_percentage){
  fixed_percentage <- as.numeric(opt$perc)
} else {
  perc_map_f <- opt$perc_map # e.g. "./map_final_perc.csv",
  percentage_map <- read.csv(perc_map_f, row.names = 1)
}

# leo matriz abundancias
list[exp, tax] <- get_abundance_and_tax_from_table(abuntable,
                                                   skip=skip_lines_abuntable)
# creo output dir
system(paste("mkdir -p", outdir))

# ===========
# creo counts
# ===========
# Solo necesito la "muestra original", para ejecutar el simulate_timeseries
if (is.null(subset_otus)) {
  counts <- exp[s]
} else {
  subset_otus <- strsplit(subset_otus, ";")[[1]]
  counts <- exp[subset_otus, s, drop=F]
}


# ==========
# simulation
# ==========
## results will be stored here:
final_abund <- list()

# relative expected abundances of each PCG.
if (fix_percentage) {
  perc <- fixed_percentage
} else {
  # picked percentages from random real replicate
  perc <- sample_n(percentage_map[percentage_map$ORIG==s, ], 1)
}

# initial total abundance
total_counts <- sum(exp[s])
# final total abundance
abun_total   <- round(total_counts * perc)

# check first if there's anything to simulate
if (total_counts == 0) {
  message(paste0("There are no detectable OTUs in initial sample ",
                 s, ". Moving to next PCG..."))
} else if (total_counts * dilution < 3) {
  stop("EXIT: 3 or less bugs will be left after diluting! Consider changing your dilution factor.")

# start simulations if everything's OK
} else {
  abund_temp <- mclapply(X = 1:no_of_simulations,
                             FUN = function(iter) {

                               # 1) simulation

                               if (abun_total == 0) { # we can't simulate anything if it's 0
                                 start <- as.count(my_transpose(counts))
                                 start <- start[order(names(start))] # this order thing is to keep indices consistent
                                 trajectory <- matrix(0, ncol=length (start), nrow = no_of_dil+1)
                                 rownames(trajectory) <- 0:no_of_dil
                                 colnames(trajectory) <- names(start)
                                 trajectory["0",]=start
                               } else {
                                 trajectory <- simulate_timeseries(counts,
                                                                   dilution = dilution,
                                                                   no_of_dil = no_of_dil,
                                                                   fixation_at = fixation_at,
                                                                   grow_step = grow_step,
                                                                   is_grow_step_a_perc = is_grow_step_a_perc,
                                                                   abun_total = round(
                                                                     total_counts *
                                                                       perc),
                                                                   keep_all_timesteps = save_all)
                               }

                               print(paste("Simulation", iter, "finished for", s))

                               return(trajectory)

                             }, mc.cores = cores)

  # save data
  message(paste0("Saving data for sample ", s, "..."))
  if (save_all == T) {
    for (timepoint in 1:(no_of_dil + 1)) {
      # pick all rows number "timepoint" from all the lists in all_abund
      # "temp" refers to a single timepoint
      temp <- lapply(abund_temp, FUN = function(traj) {
        return(traj[timepoint, , drop = F] %>% as.data.frame)}) %>%
        bind_rows %>% as_tibble # one file per timepoint

      write.csv(temp,
                file = paste0(outdir, "/simul_", outputname, "_t_", timepoint - 1, ".csv"))
    }
  } else {
    final_abund <- as.data.frame(bind_rows(abund_temp))
    write.csv(final_abund,
              file = paste0(outdir,"/simul_", outputname, ".csv"))
  }
}

end_time <- Sys.time()
print(end_time - start_time)
