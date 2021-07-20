# ---  base data IO

# packages
library(hdf5r)
library(Matrix)

#' Data frame to h5
#'
#' Data frame is converted to and saved the h5 file
#' @param df Data frame of cell annotation or gene annotation
#' @param h5 The h5 file
#' @param gr_name The group name represents the property of the data frame.
#' @importFrom hdf5r H5File h5attr "h5attr<-" h5attr_names
#' @export
#'
df_to_h5 <- function(df, h5, gr_name=NULL){
  h5df <- h5$create_group(gr_name)
  h5df[['index']] = rownames(df)
  if(ncol(df)>0){
    h5df[['colnames']] = colnames(df)
  }
  # factor to levels,charactor to levels,logical to levels
  for(k in names(df)){
    if(is.factor(df[[k]])){
      e0 <- as.integer(df[[k]]) - 1L
      e0[is.na(e0)] <- -1 # -1 is na
      h5df[[k]]<- e0 # for 0 begin

      h5df[[paste0(k,'_levels')]]<- levels(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'category'
    }
    if(is.character(df[[k]])){
      str_to_lvl <- factor(df[[k]])
      e0 <- as.integer(str_to_lvl) - 1L
      e0[is.na(e0)] <- -1 # -1 is na
      h5df[[k]]<- e0
      h5df[[paste0(k,'_levels')]]<- levels(str_to_lvl)
      h5attr(h5df[[k]], 'origin_dtype') = 'string'
    }
    if(is.logical(df[[k]])){
      h5df[[k]] <- as.integer(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'bool'
    }
    if(any(is.numeric(df[[k]]),is.integer(df[[k]]))){
      h5df[[k]] <- df[[k]]
      h5attr(h5df[[k]], 'origin_dtype') = 'number'
    }
  }
}

#' H5 to dataframe
#'
#' H5 gruop is converted to usable data frame including the observes annotation (cell annotation) and variables annotation (gene annotation)
#' @param h5df The hdf5 group saved the data frame using list mode.
#' @importFrom hdf5r H5File h5attr "h5attr<-" h5attr_names
#' @return A data frame(annotation)
#' @export
#'
h5_to_df <- function(h5df){
  df_list <- list()
  df_list[['index']] <- h5df[['index']][]
  for(k in names(h5df)){
    if(length(h5attr_names(h5df[[k]]))>0){
      df_dtype <- h5attr(h5df[[k]], 'origin_dtype')
      if(df_dtype == 'category'){
        e0 <- h5df[[k]][]
        e0[e0==-1] <- NA # -1 is na
        e0 <- e0 + 1L
        lvl <- h5df[[paste0(k,'_levels')]][]
        df_list[[k]] <- structure(.Data = e0, .Label = lvl, class = 'factor')
      }
      if(df_dtype == 'string'){
        e0 <- h5df[[k]][]
        e0[e0==-1] <- NA # -1 is na
        e0 <- e0 + 1L
        lvl <- h5df[[paste0(k,'_levels')]][]
        df_list[[k]] <- as.character(structure(.Data = e0, .Label = lvl, class = 'factor'))
      }
      if(df_dtype == 'bool'){
        df_list[[k]] <- as.logical(h5df[[k]][])
      }
      if(df_dtype == 'number'){
        df_list[[k]] <- h5df[[k]][]
      }
    }
  }
  df = as.data.frame(df_list, row.names = 'index', optional = TRUE)
  if('colnames' %in% names(h5df)){
    cnames <- h5df[['colnames']][]
    if(length(cnames) >1){
      df = df[,cnames]
    }else{
      df = df
    }
  }else{
    df=df
  }
  return(df)
}


#' Matrix to H5 format
#'
#' The matrix including the dense matrix and sparse matrix is converted to the matrix in h5 format or is stored into the h5 file.
#' @param mat The matrix object including matrix(R) and sparse matrix(Matrix package)
#' @param h5 The H5 file name that we write in
#' @param gr_name The h5 gorup store the matrix (dense matrix or sparse matrix)
#' @param save.obs.name The rownames isn't be saved(FALSE by defualt)
#' @param save.var.name The colnames isn't be saved(FALSE by defualt)
#' @importFrom hdf5r H5File h5attr
#' @importFrom methods slot
#' @export
#'
matrix_to_h5 <- function(mat, h5, gr_name = NULL, save.obs.name = FALSE, save.var.name = FALSE){
  if(!gr_name %in% names(h5)){
    h5mat = h5$create_group(gr_name)
  }
  else{
    h5mat = h5[[gr_name]]
  }
  if('dgCMatrix' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    if(save.obs.name & save.var.name){
      h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
      h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
    }
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
  }
  else if('matrix' %in% class(mat)){
    h5mat[['matrix']] <- mat
    h5mat[['dims']] <- rev(dim(mat))
    if(save.obs.name & save.var.name){
      h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
      h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
    }
    h5attr(h5mat, 'datatype') <- 'Array'
    warning(paste0(substitute(gr_name), ' is dense matrix'))#2
  }
  else if('Graph' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
  }else{
    stop('The matrix type is wrong')
  }
}

#' H5 to Matrix format
#'
#' H5 group is converted to usable matrice including dense matrices and sparse matrices.
#' Dense matrices(by R) and sparse matrices(constructed by Matrix package)
#' @param h5 The name of the group the stores the matrix data in h5 file
#' @param obs_names The observes names, such as cell names
#' @param var_names The variables names, such as gene names
#' @return dense matrix or sparse matrix
#' @importFrom hdf5r H5File h5attr
#' @export
#'
h5_to_matrix <-  function(h5mat, obs.name=NULL, var.name=NULL){
  if(all(c('obs_names','var_names') %in% names(h5mat))){
    obs.name = h5mat[['obs_names']][]
    var.name = h5mat[['var_names']][]
  }
  if(h5attr(h5mat, 'datatype') == 'SparseMatrix'){
    mat <- Matrix::sparseMatrix(i = h5mat[['indices']][],
                                p = h5mat[['indptr']][],
                                x = h5mat[['values']][],
                                dims = rev(h5mat[['dims']][]),
                                index1 = FALSE)
    if(!is.null(obs.name) & !is.null(var.name)){
      dimnames(mat) <- list(var.name, obs.name)
    }
    else{
      warning('There are no dimnames in the sparse matrix')
    }
  }
  else if(h5attr(h5mat, 'datatype') == 'Array'){
    mat <- h5mat[['matrix']][,]
    if(!is.null(obs.name) & !is.null(var.name)){
      dimnames(mat) <- list(var.name, obs.name)
    }
    else{
      warning('There are no dimnames in the matrix')
    }
  }
  return(mat)
}

# h5[['obs']] file transform obs information
to_obs_ <- function(h5){
  to_obs <- h5_to_df(h5df = h5[['obs']])
  return(to_obs)
}

# h5[['var]] file transform var information
to_var_ <- function(h5, use.raw = FALSE){
  to_var <- list()
  var = h5[['var']]
  for(v in names(var)){
    to_var[[v]] <- h5_to_df(h5df = var[[v]])
  }
  return(to_var)
}

# h5[['graphs']] file transform graphs information
to_graphs_ <- function(h5){
  to_graphs <- list()
  graphs <- h5[['graphs']]
  graphs_name <- list(knn = "RNA_nn", snn = "RNA_snn")
  for(g in names(graphs)){
    to_graphs[[graphs_name[[g]]]] = h5_to_matrix(h5mat = graphs[[g]])
  }
  return(to_graphs)
}

# h5[['data']] file transform data information
to_data_ <- function(h5){
  to_data <- list()
  data <- h5[['data']]
  for(d in names(data)){
    to_data[[d]] <- h5_to_matrix(h5mat = data[[d]])
  }
  return(to_data)
}

# h5[['layers']] file transform the layer information
to_layers_ <- function(h5){
  to_layers <- list()
  layers <- h5[['layers']]
  for(l in names(layers)){
    to_layers[[l]] <- h5_to_matrix(h5mat = layers[[l]])
  }
  return(to_layers)
}

# h5[['uns']] file transform the uns information
to_uns_ <- function(h5){
  colors_list = list()
  for(colr in grep('colors', names(h5[['uns']]), value = TRUE)){
    colors_list[[colr]] <- h5[['uns']][[colr]][]
  }
  return(colors_list)
}

# h5[['dimR']] file transform the dimr information
to_dimr_ <- function(h5){
  to_dimr <- list()
  dimr <- h5[['dimR']]
  for(DR in names(dimr)){
    to_dimr[[tolower(DR)]] <- t(dimr[[DR]][,])
    colnames(to_dimr[[tolower(DR)]]) <-  paste0(DR, '_',1:ncol(to_dimr[[tolower(DR)]]))
  }
  return(to_dimr)
}

# h5[['dimR']] file transform the dimr information of seurat object
seurat.to_dimr_ <- function(h5){
  to_dimr <- list()
  dimr <- h5[['dimR']]
  for(DR in names(dimr)){
    to_dimr[[tolower(DR)]] <- Seurat::CreateDimReducObject(embeddings = t(dimr[[DR]][,]), key = paste0(DR, "_"), assay = "RNA")
  }
  return(to_dimr)
}

# h5[['spatial']] file transforma the spatial information of seurat object
seurat.to_spatial_ <- function(h5){
  to_spatial <- h5_to_spatial(h5[['spatial']])
  return(to_spatial)
}


#--- read h5 file

#' H5 to scRNAs-seq analysis object
#'
#' Read h5 and converted h5 to the scRNA-seq analysis object
#' @param target.object Denotes which object to load. Defualt is "seurat". Available options are:
#' \itemize{
#'   \item "seurat": converted the h5 file to "seurat object".
#'   \item "singlecellexperiment": converted the h5 file to "singlecellexperiment object".
#'   \item "monocle": converted the h5 file to "monocle3 object".
#' }
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @return The single cell analysis mainstream software
#' @export
#'
read_h5 <- function(file, target.object = 'seurat', assay.name = 'RNA'){
  if(target.object == 'seurat'){
    data <- seurat_read_h5(file = file, assay.name = assay.name)
  }
  else if(target.object == 'singlecellexperiment'){
    data <- sce_read_h5(file = file, assay.name = assay.name)
  }
  else{
    stop('The output object must be specified')
  }
  return(data)
}


#--- write h5

#' The scRNAs-seq analysis objec to H5
#'
#' Write h5 and  the scRNA-seq analysis object converted to h5
#' @param object.type Denotes which object to save.Available options are:
#' \itemize{
#'   \item "seurat": converted the "seurat object" to the h5 file.
#'   \item "singlecellexperiment": converted the "singlecellexperiment object" to the h5 file.
#'   \item "monocle": converted the "monocle3 object" to the h5 file.
#' }
#' @param data The scRNA-seq analysis object data.
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @param save.graphs Default is False , determing whether to save the graph(cell-cell similarity network).
#'   seurat graph is different from scanpy graph. Their relationship are set {"distances": "knn", "connectivities": "snn"} roughly.
#' @export
#'
write_h5 <- function(data, file, object.type = 'seurat', assay.name = 'RNA', save.graphs = FALSE, save.scale=FALSE){
  if(object.type == 'seurat'){
    seurat_write_h5(seurat = data, file = file, assay.name = assay.name, save.graphs = save.graphs, save.scale = save.scale)
  }
  if(object.type == 'singlecellexperiment'){
    sce_write_h5(sce = data, file = file, assay.name = assay.name)
  }
}


#--- read_h5 part

#' The read the some data in the h5 fiel
#' @param h5 the h5 file
#' @param groups the h5 groups
#'
#' @export
read_h5part <- function(h5, groups){
  h5_list <- list()
  for(h in groups){
    h5_list[[h]] <- switch(h, data = to_data_(h5),
                               obs = to_obs_(h5),
                               var = to_var_(h5),
                               dimR = seurat.to_dimr_(h5),
                               layers = to_layers_(h5),
                               graphs = to_graphs_(h5),
                               spatial = seurat.to_spatial_(h5),
                               uns = to_uns_(h5))
  }
  return(h5_list)
}















