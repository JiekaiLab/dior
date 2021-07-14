# --- SingleCellExperiment(sce) H5 IO

# library
library(hdf5r)
library(Matrix)
library(SingleCellExperiment)

# --- SingleCellExperiment write h5 file

#' The singlecellexperiment is converted to h5 file
#'
#' Singlecellexperiment object is converted to h5 file.
#' @param sce The singlecellexperiment object.
#' @param file The h5 file.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#'
sce_write_h5 <- function(sce = NULL, file = NULL, assay.name = NULL){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(class(sce) != 'SingleCellExperiment'){
    stop('object ', substitute(sce), ' class is not SingleCellExperiment object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    sce_to_h5(sce = sce, h5 = h5, assay.name = assay.name)
    h5attr(h5, 'assay_name') <- assay.name
  },
  error = function(e){
    print(e)
  },
  finally = {
    h5$close_all()
  })
}

#' The singlecellexperiment is converted to h5 file
#'
#' Singlecellexperiment object is converted to h5 file.
#' @param sce The singlecellexperiment object.
#' @param h5 The h5 file in R.
#' @param assay.name The 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#' @importFrom SummarizedExperiment assayNames assay colData rowData "assayNames<-"
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
sce_to_h5 <- function(sce, h5, assay.name){
  h5attr(h5, 'assay_name') <- assay.name
  data <- h5$create_group('data')
  var <- h5$create_group('var')
  #--- save the matrix
  if(length(assayNames(sce))> 1){
    layers <- h5$create_group('layers')
  }
  if(!'X' %in% assayNames(sce)){
    assayNames(sce)[1] <- 'X'
    print("The first 'assayNames' defaults to 'X'")
  }
  for(as in assayNames(sce)){
    if('X' == as){
      matrix_to_h5(mat = assay(sce, as), h5 = data, gr_name = as)
    }else{
      matrix_to_h5(mat = assay(sce, as), h5 = layers, gr_name = as)
    }
  }
  #--- save cell annotation
  df_to_h5(df = colData(sce), h5 = h5, gr_name = 'obs')
  #--- save var
  df_to_h5(df = rowData(sce), h5 = var, gr_name = 'X')
  #--- save reduction dimension
  if(length(reducedDimNames(sce))>0){
    dimReduction <- h5$create_group('dimR')
    for(d in reducedDimNames(sce)){
      D = toupper(d)
      dimReduction[[D]] <- t(reducedDim(sce, d))
    }
  }

}

# --- singlecellexperiment read the h5 file

#' H5 to singlecellexperiment object
#'
#' Read h5 and converted h5 to the singlecellexperiment object
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#'
sce_read_h5 <- function(file=NULL, assay.name = NULL){
  if(!file.exists(file)){
    stop('No such file or directory')
  }
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    data <- h5_to_sce(h5 = h5, assay.name = assay.name)
  },
  error = function(e) {
    print(e)
  },
  finally = {
    h5$close_all()
  }
  )
  return(data)
}



#' H5 to the singlecellexperiment
#'
#' H5 file is converted to the singlecellexperiment obejct
#' @param h5 The h5 file in R
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#'
h5_to_sce <- function(h5, assay.name = 'RNA',use.data ='X'){
  options(warn = -1)
  if(h5attr(h5, 'assay_name') == assay.name){
    sce_list <- list()
    for(h in names(h5)){
      # tran to the list
      sce_list[[h]] <- switch(h, data = to_data_(h5),
                              obs = to_obs_(h5),
                              var = to_var_(h5),
                              dimR = to_dimr_(h5),
                              layers = to_layers_(h5),
                              graphs = to_graphs_(h5),
                              uns = to_uns_(h5))
    }
    if('dimR' %in% names(sce_list)){
      dimR_list <- sce_list[['dimR']]
      for(dr in names(dimR_list)){
        rownames(sce_list[['dimR']][[dr]]) <- rownames(sce_list[['obs']])
      }
    }
    if(use.data == 'rawX'){
      ass <- sce_list[['data']][use.data]
      rownames(ass[[use.data]]) <-  rownames(sce_list[['var']][[use.data]])
      colnames(ass[[use.data]]) <- rownames(sce_list[['obs']])
    }
    if(use.data == 'X'){
      ass <- list(X = sce_list[['data']][[use.data]])
      rownames(ass[[use.data]]) <-  rownames(sce_list[['var']][[use.data]])
      colnames(ass[[use.data]]) <- rownames(sce_list[['obs']])
      if('layers' %in% names(sce_list)){
        for(l in names(sce_list[['layers']])){
          ass[[l]] <- sce_list[['layers']][[l]]
        }
        for(d in names(ass)){
          rownames(ass[[d]]) <-  rownames(sce_list[['var']][[use.data]])
          colnames(ass[[d]]) <- rownames(sce_list[['obs']])
        }
      }
    }
    # save the all dataset and doing
    sce <- SingleCellExperiment::SingleCellExperiment(assays =ass,
                                                      colData = sce_list[['obs']],
                                                      rowData = sce_list[['var']][[use.data]],
                                                      reducedDims = sce_list[['dimR']])

  }
  return(sce)
}






