# group by PCG
# source: https://stackoverflow.com/questions/45720062/grouping-igraph-vertices-in-a-weighted-network-by-color-subgroup-in-r/45802982#45802982?newreg=81fceb02f64f42dbadc14fded775b529
GroupByPCG1 = function(Groups, spacing = 5) {
  Position = (order(Groups) + spacing*Groups)
  Angle    = Position * 2 * pi / max(Position)
  matrix(c(cos(Angle), sin(Angle)), ncol=2)
}

GroupByPCG2 = function(Groups, spacing = 5) { 
  numGroups = length(unique(Groups))
  GAngle    = (1:numGroups) * 2 * pi / numGroups
  Centers   = matrix(c(cos(GAngle), sin(GAngle)), ncol=2)
  x = y = c()
  for(i in 1:numGroups) {
    curGroup = which(Groups == unique(Groups)[i])
    VAngle = (1:length(curGroup)) * 2 * pi / length(curGroup)
    x = c(x, Centers[i,1] + cos(VAngle) / numGroups )
    y = c(y, Centers[i,2] + sin(VAngle) / numGroups)
  }
  matrix(c(x, y), ncol=2)
}

my_transpose <- function(df){
  library(data.table,quietly = 1)
  t_df=data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  return(t_df)
}

#https://stackoverflow.com/questions/10689055/create-an-empty-data-frame
create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}


my_plot_lines <-  function(a, colors, limit=0,ylim_=0.15) {
  if (limit!=0) {j=limit} else {j=dim(a)[2]}
  
  plot(t(a)[1:j,1],type = "l", ylim = c(0.00,ylim_),col=colors[1])
  
  # for (i in c(2,10, #13 # TODO remove this
  #             12,#16
  #             6,3)){
  
  for (i in 2:dim(a)[1]){
    lines(x=as.vector(t(a[i,1:j])),
          col=colors [i])
  }
}

get_abundance_and_tax_from_table <- function (exp_file, # required
                                              species_are_rows=TRUE,
                                              sep="\t", # extra internal options you can change
                                              row.names=1,
                                              skip=1,
                                              taxa_fields=NULL,
                                              tax_col_name="taxonomy",
                                              tax_sep=";",# 4 chars or longer "breaks" renamer (will include this as if it were a real name)
                                              NA_option="___",
                                              check.names=FALSE) { 
  # This function reads a .csv/.tsv speciesXsamples file, where the last "sample" 
  # column (or row) is the taxonomy. Then returns the abundances and taxa data in
  # separate data.frames
  require(gsubfn)
  require(tidyverse)
  if (is.null(taxa_fields)) {taxa_fields = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")}
  
  exp=read.csv(exp_file,sep = sep,skip = skip,row.names=row.names,check.names = check.names)
  if (!species_are_rows) {exp<-my_transpose(exp)}
  
  tax<-exp["taxonomy"]; exp<-exp[1:dim(exp)[2]-1]
  tax <- tax %>% separate("taxonomy",sep = tax_sep,taxa_fields)
  tax[is.na(tax)]<- NA_option # avoid na-related errors
  
  return(list(exp,tax))
}

filter_ASVs_by_presence_in_perc_samples <- function(ac_full,perc=0.10,ASVs_are_rows=TRUE){ # TODO include junk,
  # me quedo solo con las presentes en el 10% o más de las muestras
  if (ASVs_are_rows==FALSE) {
    ac_full <- my_transpose(ac_full)
  }
  
  total_samples <-dim(ac_full)[2]
  acpresence <- ac_full
  acpresence [acpresence!=0] <-1
  ac1      <- ac_full [rowSums(acpresence)>=perc*total_samples,]
  
  if (ASVs_are_rows==FALSE) {
    ac1 <- my_transpose(ac1)
  }
  return(ac1)
}


filter_ASVs_by_relative_abundance_in_perc_samples <- function(ac_full,rel_abun=0.01,perc=0.10,ASVs_are_rows=TRUE, include_junk=FALSE){
  # me quedo con las que tengan al menos una abundancia relativa (al total) del 1% en al menos 10%
  # (tienen al menos 200 reads en al menos 9 muestras)
  if (ASVs_are_rows==FALSE) {
    ac_full <- my_transpose(ac_full)
  }
  
  total_samples <-dim(ac_full)[2]
  total_per_sample<-colSums(ac_full)[1]
  acabuncheck <- ac_full>=rel_abun*total_per_sample # para cada muestra para cada read, llega a 200?
  ac2   <- ac_full [rowSums(acabuncheck)>=perc*total_samples,]
  
  if (include_junk){
    ac2=rbind(ac2,(total_per_sample-colSums(ac2))) # Para que el total por sample siga valiendo lo mismo
    rownames(ac2)[dim(ac2)[1]] <- "junk"
  }
  
  if (ASVs_are_rows==FALSE) {
    ac2 <- my_transpose(ac2)
  }
  return(ac2)
}


# esta la creo para quedarme con aquellas que aparezcan en más de 2 samples 
# de mi muestra inicial, haya 6,7 u 8 en total.
filter_ASVs_by_relative_abundance_in_num_samples <- function(ac_full,rel_abun=0.01,num=2,ASVs_are_rows=TRUE, include_junk=FALSE){
  # me quedo con las que tengan al menos una abundancia relativa (al total) del 1% en al menos 2 muestras
  if (ASVs_are_rows==FALSE) {
    ac_full <- my_transpose(ac_full)
  }
  
  total_samples <-dim(ac_full)[2]
  total_per_sample<-colSums(ac_full)[1]
  acabuncheck <- ac_full>=rel_abun*total_per_sample # para cada muestra para cada read, llega a 200?
  ac2   <- ac_full [rowSums(acabuncheck)>=num,]
  
  if (include_junk){
    ac2=rbind(ac2,(total_per_sample-colSums(ac2))) # Para que el total por sample siga valiendo lo mismo
    rownames(ac2)[dim(ac2)[1]] <- "junk"
  }
  
  if (ASVs_are_rows==FALSE) {
    ac2 <- my_transpose(ac2)
  }
  return(ac2)
}


filter_taxa_by_relative_abundance_in_perc_samples <- function(tax_full,ac_full,rel_abun=0.01,perc=0.10,taxa_are_rows=TRUE,ASVs_are_rows=TRUE){
  # me quedo con las que tengan al menos una abundancia relativa (al total) del 1% en al menos 10%
  # (tienen al menos 200 reads en al menos 9 muestras)
  if (taxa_are_rows==FALSE) {
    tax_full <- my_transpose(tax_full)
  }
  if (ASVs_are_rows==FALSE) {
    ac_full <- my_transpose(ac_full)
  }
  total_samples <-dim(ac_full)[2]
  total_per_sample<-colSums(ac_full)[1]
  acabuncheck <- ac_full>=rel_abun*total_per_sample # para cada muestra para cada read, llega a 200?
  tax2   <- tax_full [rowSums(acabuncheck)>=perc*total_samples,]
  
  if (taxa_are_rows==FALSE) {
    tax2 <- my_transpose(tax2)
  }
  return(tax2)
}



filter_taxa_by_relative_abundance_in_num_samples <- function(tax_full,ac_full,rel_abun=0.01,num=2,taxa_are_rows=TRUE,ASVs_are_rows=TRUE){
  # me quedo con las que tengan al menos una abundancia relativa (al total) del 1% en al menos 2 muestras
  if (taxa_are_rows==FALSE) {
    tax_full <- my_transpose(tax_full)
  }
  if (ASVs_are_rows==FALSE) {
    ac_full <- my_transpose(ac_full)
  }
  total_samples <-dim(ac_full)[2]
  total_per_sample<-colSums(ac_full)[1]
  acabuncheck <- ac_full>=rel_abun*total_per_sample # para cada muestra para cada read, llega a 200?
  tax2   <- tax_full [rowSums(acabuncheck)>=num,]
  
  if (taxa_are_rows==FALSE) {
    tax2 <- my_transpose(tax2)
  }
  return(tax2)
}

## variables
colorsmaker <-function() {
  colorsmaker = rep(c('violetred', 'orchid', "tomato", "thistle", "slategray", 
                      "deepskyblue", "gold", "blue",  "darkolivegreen", "yellow", "springgreen",
                      "goldenrod", "orange", "steelblue", "lightcoral",
                      "darkred", "brown", "white"),7)
}

colorsmaker2 <-function() {
  colorsmaker = rep(c("darkred", "brown", "firebrick", "crimson", "red", "tomato", "coral", "indianred", "lightcoral", "darksalmon", "salmon", "lightsalmon", "orangered", "darkorange", "orange", "gold", "darkgoldenrod", "goldenrod", "palegoldenrod", "darkkhaki", "khaki", "olive", "yellow", "yellowgreen", "darkolivegreen", "olivedrab", "lawngreen", "chartreuse", "greenyellow", "darkgreen", "green", "forestgreen", "lime", "limegreen", "lightgreen", "palegreen", "darkseagreen", "mediumspringgreen", "springgreen", "seagreen", "mediumaquamarine", "mediumseagreen", "lightseagreen", "darkslategray", "teal", "darkcyan", "aqua", "cyan", "lightcyan", "darkturquoise", "turquoise", "mediumturquoise", "paleturquoise", "aquamarine", "powderblue", "cadetblue", "steelblue", "cornflowerblue", "deepskyblue", "dodgerblue", "lightblue", "skyblue", "lightskyblue", "midnightblue", "navy", "darkblue", "mediumblue", "blue", "royalblue", "blueviolet", "indigo", "darkslateblue", "slateblue", "mediumslateblue", "mediumpurple", "darkmagenta", "darkviolet", "darkorchid", "mediumorchid", "purple", "thistle", "plum", "violet", "magenta", "orchid", "mediumvioletred", "palevioletred", "deeppink", "hotpink", "lightpink", "pink", "antiquewhite", "beige", "bisque", "blanchedalmond", "wheat", "cornsilk", "lemonchiffon", "lightgoldenrodyellow", "lightyellow", "saddlebrown", "sienna", "chocolate", "peru", "sandybrown", "burlywood", "tan", "rosybrown", "moccasin", "navajowhite", "peachpuff", "mistyrose", "lavenderblush", "linen", "oldlace", "papayawhip", "seashell", "mintcream", "slategray", "lightslategray", "lightsteelblue", "lavender", "floralwhite", "aliceblue", "ghostwhite", "honeydew", "ivory", "azure", "snow", "black", "dimgrey", "grey", "darkgrey", "silver", "lightgrey", "gainsboro", "whitesmoke", "white"),7)
}

colorsmaker_hex <- function(){
  colorsmaker_hex = c("#B0171F", "#DC143C", "#FFB6C1", "#FFAEB9", "#EEA2AD", "#CD8C95", "#8B5F65", "#FFC0CB", "#FFB5C5", "#EEA9B8", "#CD919E", "#8B636C", "#DB7093", "#FF82AB", "#EE799F", "#CD6889", "#8B475D", "#FFF0F5", "#EEE0E5", "#CDC1C5", "#8B8386", "#FF3E96", "#EE3A8C", "#CD3278", "#8B2252", "#FF69B4", "#FF6EB4", "#EE6AA7", "#CD6090", "#8B3A62", "#872657", "#FF1493", "#EE1289", "#CD1076", "#8B0A50", "#FF34B3", "#EE30A7", "#CD2990", "#8B1C62", "#C71585", "#D02090", "#DA70D6", "#FF83FA", "#EE7AE9", "#CD69C9", "#8B4789", "#D8BFD8", "#FFE1FF", "#EED2EE", "#CDB5CD", "#8B7B8B", "#FFBBFF", "#EEAEEE", "#CD96CD", "#8B668B", "#DDA0DD", "#EE82EE", "#FF00FF", "#EE00EE", "#CD00CD", "#8B008B", "#800080", "#BA55D3", "#E066FF", "#D15FEE", "#B452CD", "#7A378B", "#9400D3", "#9932CC", "#BF3EFF", "#B23AEE", "#9A32CD", "#68228B", "#4B0082", "#8A2BE2", "#9B30FF", "#912CEE", "#7D26CD", "#551A8B", "#9370DB", "#AB82FF", "#9F79EE", "#8968CD", "#5D478B", "#483D8B", "#8470FF", "#7B68EE", "#6A5ACD", "#836FFF", "#7A67EE", "#6959CD", "#473C8B", "#F8F8FF", "#E6E6FA", "#0000FF", "#0000EE", "#0000CD", "#00008B", "#000080", "#191970", "#3D59AB", "#4169E1", "#4876FF", "#436EEE", "#3A5FCD", "#27408B", "#6495ED", "#B0C4DE", "#CAE1FF", "#BCD2EE", "#A2B5CD", "#6E7B8B", "#778899", "#708090", "#C6E2FF", "#B9D3EE", "#9FB6CD", "#6C7B8B", "#1E90FF", "#1C86EE", "#1874CD", "#104E8B", "#F0F8FF", "#4682B4", "#63B8FF", "#5CACEE", "#4F94CD", "#36648B", "#87CEFA", "#B0E2FF", "#A4D3EE", "#8DB6CD", "#607B8B", "#87CEFF", "#7EC0EE", "#6CA6CD", "#4A708B", "#87CEEB", "#00BFFF", "#00B2EE", "#009ACD", "#00688B", "#33A1C9", "#ADD8E6", "#BFEFFF", "#B2DFEE", "#9AC0CD", "#68838B", "#B0E0E6", "#98F5FF", "#8EE5EE", "#7AC5CD", "#53868B", "#00F5FF", "#00E5EE", "#00C5CD", "#00868B", "#5F9EA0", "#00CED1", "#F0FFFF", "#E0EEEE", "#C1CDCD", "#838B8B", "#E0FFFF", "#D1EEEE", "#B4CDCD", "#7A8B8B", "#BBFFFF", "#AEEEEE", "#96CDCD", "#668B8B", "#2F4F4F", "#97FFFF", "#8DEEEE", "#79CDCD", "#528B8B", "#00FFFF", "#00EEEE", "#00CDCD", "#008B8B", "#008080", "#48D1CC", "#20B2AA", "#03A89E", "#40E0D0", "#808A87", "#00C78C", "#7FFFD4", "#76EEC6", "#66CDAA", "#458B74", "#00FA9A", "#F5FFFA", "#00FF7F", "#00EE76", "#00CD66", "#008B45", "#3CB371", "#54FF9F", "#4EEE94", "#43CD80", "#2E8B57", "#00C957", "#BDFCC9", "#3D9140", "#F0FFF0", "#E0EEE0", "#C1CDC1", "#838B83", "#8FBC8F", "#C1FFC1", "#B4EEB4", "#9BCD9B", "#698B69", "#98FB98", "#9AFF9A", "#90EE90", "#7CCD7C", "#548B54", "#32CD32", "#228B22", "#00FF00", "#00EE00", "#00CD00", "#008B00", "#008000", "#006400", "#308014", "#7CFC00", "#7FFF00", "#76EE00", "#66CD00", "#458B00", "#ADFF2F", "#CAFF70", "#BCEE68", "#A2CD5A", "#6E8B3D", "#556B2F", "#6B8E23", "#C0FF3E", "#B3EE3A", "#9ACD32", "#698B22", "#FFFFF0", "#EEEEE0", "#CDCDC1", "#8B8B83", "#F5F5DC", "#FFFFE0", "#EEEED1", "#CDCDB4", "#8B8B7A", "#FAFAD2", "#FFFF00", "#EEEE00", "#CDCD00", "#8B8B00", "#808069", "#808000", "#BDB76B", "#FFF68F", "#EEE685", "#CDC673", "#8B864E", "#F0E68C", "#EEE8AA", "#FFFACD", "#EEE9BF", "#CDC9A5", "#8B8970", "#FFEC8B", "#EEDC82", "#CDBE70", "#8B814C", "#E3CF57", "#FFD700", "#EEC900", "#CDAD00", "#8B7500", "#FFF8DC", "#EEE8CD", "#CDC8B1", "#8B8878", "#DAA520", "#FFC125", "#EEB422", "#CD9B1D", "#8B6914", "#B8860B", "#FFB90F", "#EEAD0E", "#CD950C", "#8B6508", "#FFA500", "#EE9A00", "#CD8500", "#8B5A00", "#FFFAF0", "#FDF5E6", "#F5DEB3", "#FFE7BA", "#EED8AE", "#CDBA96", "#8B7E66", "#FFE4B5", "#FFEFD5", "#FFEBCD", "#FFDEAD", "#EECFA1", "#CDB38B", "#8B795E", "#FCE6C9", "#D2B48C", "#9C661F", "#FF9912", "#FAEBD7", "#FFEFDB", "#EEDFCC", "#CDC0B0", "#8B8378", "#DEB887", "#FFD39B", "#EEC591", "#CDAA7D", "#8B7355", "#FFE4C4", "#EED5B7", "#CDB79E", "#8B7D6B", "#E3A869", "#ED9121", "#FF8C00", "#FF7F00", "#EE7600", "#CD6600", "#8B4500", "#FF8000", "#FFA54F", "#EE9A49", "#CD853F", "#8B5A2B", "#FAF0E6", "#FFDAB9", "#EECBAD", "#CDAF95", "#8B7765", "#FFF5EE", "#EEE5DE", "#CDC5BF", "#8B8682", "#F4A460", "#C76114", "#D2691E", "#FF7F24", "#EE7621", "#CD661D", "#8B4513", "#292421", "#FF7D40", "#FF6103", "#8A360F", "#A0522D", "#FF8247", "#EE7942", "#CD6839", "#8B4726", "#FFA07A", "#EE9572", "#CD8162", "#8B5742", "#FF7F50", "#FF4500", "#EE4000", "#CD3700", "#8B2500", "#5E2612", "#E9967A", "#FF8C69", "#EE8262", "#CD7054", "#8B4C39", "#FF7256", "#EE6A50", "#CD5B45", "#8B3E2F", "#8A3324", "#FF6347", "#EE5C42", "#CD4F39", "#8B3626", "#FA8072", "#FFE4E1", "#EED5D2", "#CDB7B5", "#8B7D7B", "#FFFAFA", "#EEE9E9", "#CDC9C9", "#8B8989", "#BC8F8F", "#FFC1C1", "#EEB4B4", "#CD9B9B", "#8B6969", "#F08080", "#CD5C5C", "#FF6A6A", "#EE6363", "#8B3A3A", "#CD5555", "#A52A2A", "#FF4040", "#EE3B3B", "#CD3333", "#8B2323", "#B22222", "#FF3030", "#EE2C2C", "#CD2626", "#8B1A1A", "#FF0000", "#EE0000", "#CD0000", "#8B0000", "#800000", "#8E388E", "#7171C6", "#7D9EC0", "#388E8E", "#71C671", "#8E8E38", "#C5C1AA", "#C67171", "#555555", "#1E1E1E", "#282828", "#515151", "#5B5B5B", "#848484", "#8E8E8E", "#AAAAAA", "#B7B7B7", "#C1C1C1", "#EAEAEA", "#F4F4F4", "#FFFFFF", "#F5F5F5", "#DCDCDC", "#D3D3D3", "#C0C0C0", "#A9A9A9", "#808080", "#696969", "#000000", "#800000", "#8B0000", "#A52A2A", "#B22222", "#DC143C", "#FF0000", "#FF6347", "#FF7F50", "#CD5C5C", "#F08080", "#E9967A", "#FA8072", "#FFA07A", "#FF4500", "#FF8C00", "#FFA500", "#FFD700", "#B8860B", "#DAA520", "#EEE8AA", "#BDB76B", "#F0E68C", "#808000", "#FFFF00", "#9ACD32", "#556B2F", "#6B8E23", "#7CFC00", "#7FFF00", "#ADFF2F", "#006400", "#008000", "#228B22", "#00FF00", "#32CD32", "#90EE90", "#98FB98", "#8FBC8F", "#00FA9A", "#00FF7F", "#2E8B57", "#66CDAA", "#3CB371", "#20B2AA", "#2F4F4F", "#008080", "#008B8B", "#00FFFF", "#00FFFF", "#E0FFFF", "#00CED1", "#40E0D0", "#48D1CC", "#AFEEEE", "#7FFFD4", "#B0E0E6", "#5F9EA0", "#4682B4", "#6495ED", "#00BFFF", "#1E90FF", "#ADD8E6", "#87CEEB", "#87CEFA", "#191970", "#000080", "#00008B", "#0000CD", "#0000FF", "#4169E1", "#8A2BE2", "#4B0082", "#483D8B", "#6A5ACD", "#7B68EE", "#9370DB", "#8B008B", "#9400D3", "#9932CC", "#BA55D3", "#800080", "#D8BFD8", "#DDA0DD", "#EE82EE", "#FF00FF", "#DA70D6", "#C71585", "#DB7093", "#FF1493", "#FF69B4", "#FFB6C1", "#FFC0CB", "#FAEBD7", "#F5F5DC", "#FFE4C4", "#FFEBCD", "#F5DEB3", "#FFF8DC", "#FFFACD", "#FAFAD2", "#FFFFE0", "#8B4513", "#A0522D", "#D2691E", "#CD853F", "#F4A460", "#DEB887", "#D2B48C", "#BC8F8F", "#FFE4B5", "#FFDEAD", "#FFDAB9", "#FFE4E1", "#FFF0F5", "#FAF0E6", "#FDF5E6", "#FFEFD5", "#FFF5EE", "#F5FFFA", "#708090", "#778899", "#B0C4DE", "#E6E6FA", "#FFFAF0", "#F0F8FF", "#F8F8FF", "#F0FFF0", "#FFFFF0", "#F0FFFF", "#FFFAFA", "#000000", "#696969", "#808080", "#A9A9A9", "#C0C0C0", "#D3D3D3",          
                      "#DCDCDC", "#F5F5F5", "#FFFFFF")
}

pcgcodemaker <- function(tax,f_=FALSE) {
  # Entero == 1; Pseudo == 2; others == 3
  if (f_==TRUE){
    e=" f__Enterobacteriaceae"
    p=" f__Pseudomonadaceae"
  } else {
    e="Enterobacteriaceae"
    p="Pseudomonadaceae"
  }
  entero000 = as.numeric(tax["Family"]!=e)*2+1
  pseudo000 = as.numeric(tax["Family"]==p)
  pcgcode= factor(entero000-pseudo000, ordered=TRUE); rm(entero000); rm(pseudo000)
  return(pcgcode)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

pcgcodemaker2 <- function(tax,together=TRUE) {
  if (together==TRUE){
    b=c(" o__Burkholderiales", " o__Xanthomonadales")
    p=" o__Pseudomonadales"
  } else {
    b=" o__Burkholderiales"
    x=" o__Xanthomonadales"
    p=" o__Pseudomonadales"
  }
  b000 = as.numeric(sapply(tax["Order"] ,FUN=function(x){x%!in% b}))*2+1
  pseudo000 = as.numeric(tax["Order"]==p)
  pcgcode= b000-pseudo000
  
  if (together){
    pcgcode <- factor(as.numeric(pcgcode,ordered = TRUE))
  } else {
    x000 = as.numeric(tax["Order"]==x)
    pcgcode[x000==1]<-4  # x 4, b 1, p 2, others 3
    pcgcode <- factor(as.numeric(pcgcode,ordered = TRUE))
  }
  return(pcgcode)
}

add_pcgcode_legend <- function(x="topright",legend = c("entero", "pseudo", "others"),col = c(1,2,3), cex=0.9, lty = 1){
  legend(x=x,legend = legend,col=col, cex=cex, lty=lty)
}

add_colorcode_legend <- function(x="topright",tax = NULL, cex=0.9, lty = 1){
  legend = rownames(tax)
  col    = colorsmaker()[1:dim(tax)[1]]
  legend(x=x,legend = legend,col=col, cex=cex, lty=lty)
}



# Robada de https://github.com/shkonishi/cornet/blob/master/R/matoedge.R
matoedge <- function(mat, diag = FALSE, zero.weight = FALSE, format="igraph"){
  # argument check: mat is a squared matrix
  if(nrow(mat) != ncol(mat)){
    stop("This matrix is not squared matrix")
  }
  
  if (format == "igraph"){
    if (diag == FALSE){
      igraph::graph.adjacency(mat, mode="undirected", weighted=TRUE, diag = FALSE)
    } else {
      igraph::graph.adjacency(mat, mode="undirected", weighted=TRUE)
    }
    
  } else if (format =="df"){
    if (diag == FALSE){
      bit <- cbind(rep(1:(nrow(mat)-1), (nrow(mat)-1):1),
                   unlist(lapply(seq(nrow(mat))[-1], function(i)i:(nrow(mat)))) )
      
      if(identical(rownames(mat), NULL)){
        x_id=bit[,1]; y_id=bit[,2]
      } else {
        x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
      }
      d <- data.frame(x_id, y_id, value = mat[lower.tri(mat)])
      if(zero.weight==FALSE){ d[d$value != 0,] } else if (zero.weight==TRUE){d}
      
    }else{
      bit <- cbind(rep(1:nrow(mat), nrow(mat):1),
                   unlist(lapply(1:ncol(mat), function(i)seq(i, ncol(mat)))))
      
      if(identical(rownames(mat), NULL)){
        x_id=bit[,1]; y_id=bit[,2]
      } else {
        x_id=rownames(mat)[bit[,1]]; y_id=rownames(mat)[bit[,2]]
      }
      d <- data.frame(x_id, y_id, value = mat[lower.tri(mat, diag = T)])
      if(zero.weight==FALSE){ d[d$value != 0,] } else if (zero.weight==TRUE){d}
      
    }
  }
}

renamer=function(tax, width){
  library(stringr,quietly = 1)
  final_name=list()
  for (otu in rownames(tax)) {
    final_name[[otu]]=as.character(otu)
    for (tax_level in colnames(tax)[5:7]){
      taxon <- tax[otu,tax_level]
      if (!is.na(taxon)) {
        if (nchar(taxon)-3 < width) {
          newwidth <- nchar(taxon)-3
        } else {
          newwidth <- width
        }
        tr<-str_split(taxon,"__")[[1]][2]
        tr<-str_trunc(tr,ellipsis = "", width = newwidth)
        final_name[[otu]]<-paste0(final_name[[otu]],"_",tr)
        final_name[[otu]] <- gsub("_+$", "", final_name[[otu]])
      }
    }
  }
  return(final_name)
}

prettyrenamer <- function(tax, taxtoinclude = 1){
  levels <- c("(K)", "(P)", "(C)", "(O)", "(F)", "(G)", "(S)")
  final_name <- list()
  for (otu in rownames(tax)) {
    lvl <- 0
    tax_list <- c()
    for (tax_level in colnames(tax)){
      taxon <- tax[otu,tax_level]
      if (!is.na(taxon)) {
        if (nchar(taxon) > 4) {
          lvl <- lvl + 1
          tr<-str_split(taxon, "__")[[1]][2]
          tr<-str_trunc(tr, ellipsis = "", width = nchar(taxon)-3)
          tax_list <- c(tax_list, tr)
        }
      }
    }
    if (lvl==7) {tti <- taxtoinclude + 1} else {tti <- taxtoinclude}# if species, I should include the genus too
    final_name[[otu]] <- paste(tax_list[(lvl-(tti-1)):lvl], collapse = "_")
    final_name[[otu]] <- paste(as.character(otu), final_name[[otu]], levels[lvl], sep = "_")
  }
  return(final_name)
}

sa_renamer_glc = function(input, sa_as_cols = TRUE) {
  # ugly function.  probably many better ways to make it work.
  # es asumiendo que las SRR no están por orden, lo cual no va a pasar nunca.
  of  <- '/home/urihs/AAA/2021-05-14__estudio_muestras_iniciales/orig_table.from_biom_0.99.txt'
  ofr <- read.csv(of,sep="\t",row.names = 1, skip=1)
  map_of <- "map_glucosa_orig.csv"
  map_ofr<- read.csv(map_of,sep=",", row.names = 1)
  
  f= "/home/urihs/AAA/OTU_data/glucosa/0.99/table_glucosa.txt"
  fr <- read.csv(f,sep="\t",row.names = 1, skip=1)
  map_f <- "/home/silvia/AAA/data/map_glc.csv"
  map_fr<- read.csv(map_f,sep=",")
  
  # filtrar por finales
  map_ofr <- map_ofr[map_ofr$TIME=="final",]
  
  # ordenar
  map_fr    <- map_fr [order(map_fr$SRR) ,]
  map_ofr   <- map_ofr[order(map_ofr$SRR),]
  
  # crear índice
  index     <- cbind(map_ofr, map_fr$SA)
  colnames(index["map_fr$SA"])<-"new"
  
  # renombrar
  if (sa_as_cols) {
    new <-index[colnames(input),]$new
  }
  return(new)
}

# https://datascience.stackexchange.com/questions/70097/find-the-mode-value-and-frequency-in-r
mode <- function(x) {
  if ( anyNA(x) ) x = x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# newname="otu%s"
# rownames(ac_full) <-sprintf(newname,seq(1:dim(ac_full)[1]))


# # normalization 
# pseudo=0.5
# ac2.norm <- ac2 + pseudo;
# ac2.norm <-ac2.norm / colSums(ac_full)

# tax <- tax %>% separate("taxonomy",sep = ";",c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))




## Quitar el .fastq.gz
# SRRlist=c()
# for (srr in temp){
#   SRRlist<-c(SRRlist,strsplit(srr,split=".fa",fixed = F)[[1]][1])
# }






# colors
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA",
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")