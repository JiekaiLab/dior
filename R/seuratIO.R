# Seurat object H5 IO

# library
library(Seurat)
library(hdf5r)
library(Matrix)
library(Hmisc)

# --- seurat read the H5 file

# seurat insert the dimr information
seurat.insert_dimr_ <- function(seurat, seurat_list){
  for(dr in names(seurat_list[['dimR']])){
    seurat@reductions[[dr]] <- seurat_list[['dimR']][[dr]]
    rownames(seurat@reductions[[dr]]@cell.embeddings) = rownames(seurat[[]])
  }
  return(seurat)
}

# seurat insert the spatial information
seurat.insert_spatial_ <- function(seurat, seurat_list){
  seurat@images <- seurat_list[['spatial']]
  return(seurat)
}

# seurat insert the graphs information
seurat.insert_graphs_ <- function(seurat, seurat_list){
  for(g in names(seurat_list[['graphs']])){
    rownames(seurat_list[['graphs']][[g]]) <- rownames(seurat[[]])
    colnames(seurat_list[['graphs']][[g]]) <- rownames(seurat[[]])
    seurat@graphs[[g]] <- seurat_list[['graphs']][[g]]
  }
  return(seurat)
}

# seurat insert the uns information
seurat.insert_uns_ <- function(seurat, seurat_list){
  seurat@misc <- seurat_list[['uns']]
  return(seurat)
}

# seurat insert the layers information
seurat.insert_layers_ <- function(seurat, seurat_list){
  for(l in names(seurat_list[['layers']])){
    rownames(seurat_list[['layers']][[l]]) <- rownames(seurat_list[['var']][['X']])
    colnames(seurat_list[['layers']][[l]]) <- rownames(seurat[[]])
    seurat@assays[[l]] <- Seurat::CreateAssayObject(counts = seurat_list[['layers']][[l]])
  }
  return(seurat)
}

#' H5 to Seuart object
#'
#' Read h5 and converted h5 to the seurat object
#' @param file The h5 file
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#'
seurat_read_h5 <- function(file=NULL, assay.name = NULL){
  if(!file.exists(file)){
    stop('No such file or directory')
  }
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    data <- h5_to_seurat(h5 = h5, assay.name = assay.name)
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

#' H5 to the seurat
#'
#' H5 file is converted to the seurat obejct
#' @param h5 The h5 file in R
#' @param assay.name 'assay.name' is used to flag the data type, and the default is "RNA", meaning this is scRNA-seq data.
#' @importFrom hdf5r H5File h5attr
#'
h5_to_seurat <- function(h5, assay.name){
  options(warn = -1)
  if(h5attr(h5, 'assay_name') == assay.name){
    seurat_list <- list()
    for(h in names(h5)){
      seurat_list[[h]] <- switch(h, data = to_data_(h5),
                                 obs = to_obs_(h5),
                                 var = to_var_(h5),
                                 dimR = seurat.to_dimr_(h5),
                                 layers = to_layers_(h5),
                                 graphs = to_graphs_(h5),
                                 spatial = seurat.to_spatial_(h5),
                                 uns = to_uns_(h5))
    }
  }
  #--- create the seurat object: add the counts or norm data
  if(assay.name == 'spatial'){
    assay.name <- capitalize(assay.name)
  }
  rownames(seurat_list[['data']][['X']]) <- rownames(seurat_list[['var']][['X']])
  colnames(seurat_list[['data']][['X']]) <- rownames(seurat_list[['obs']])
  if(all(c('X', 'rawX') %in% names(seurat_list[['data']]))){
    rownames(seurat_list[['data']][['rawX']]) <-  rownames(seurat_list[['var']][['rawX']])
    colnames(seurat_list[['data']][['rawX']]) <- rownames(seurat_list[['obs']])
    seurat <- Seurat::CreateSeuratObject(counts = seurat_list[['data']][['rawX']], assay = assay.name)
    seurat@assays[[assay.name]]@scale.data <- seurat_list[['data']][['X']]
    seurat@assays[[assay.name]]@meta.features <- seurat_list[['var']][['rawX']]
  }else{
    seurat <- Seurat::CreateSeuratObject(counts = seurat_list[['data']][['X']], assay = assay.name)
    seurat@assays[[assay.name]]@meta.features <- seurat_list[['var']][['X']]
  }
  seurat@meta.data <- seurat_list[['obs']]
at@misc <- seurat_list[['uns']]
  for(sl in names(seurat_list)){
    seurat <- switch(sl,
                     data = seurat,
                     obs = seurat,
                     var = seurat,
                     dimR = seurat.insert_dimr_(seurat, seurat_list),
                     layers = seurat.insert_layers_(seurat, seurat_list),
                     graphs = seurat.insert_graphs_(seurat, seurat_list),
                     spatial = seurat.insert_spatial_(seurat, seurat_list),
                     uns = seurat.insert_uns_(seurat, seurat_list))
  }
  return(seurat)
}

# --- seurat write the h5 file

#' The seurat object to h5
#'
#' The Seurat object is converted to the h5 file.
#' @param seurat The seurat object.
#' @param file The h5 file.
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @param save.graphs Default is False , determing whether to save the graph(cell-cell similarity network).
#'   seurat graph is different from scanpy graph. Their relationship are set {"distances": "knn", "connectivities": "snn"} roughly.
#' @importFrom hdf5r H5File h5attr
#' @export
#'
seurat_write_h5 <- function(seurat = NULL, file = NULL, assay.name = NULL, save.graphs = FALSE, save.scale = FALSE){
  if(is.null(file)){
    stop('No such file or directory')
  }
  if(class(seurat) != 'Seurat'){
    stop('object ', substitute(seurat), ' class is not Seurat object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    seurat_to_h5(seurat = seurat, h5 = h5, assay.name = assay.name, save.graphs = save.graphs, save.scale = save.scale)
    h5attr(h5, 'assay_name') <- assay.name
  },
  error = function(e) {
    print(e)
  },
  finally = {
    h5$close_all()
  }
  )
}


#' Seurat is converted to h5 file
#'
#' Seurat object is converted to h5 file.
#' @param seurat The seurat object.
#' @param h5 The h5 file in R.
#' @param assay.name 'assay.name' is used to flag the data type. Defualt is "RNA". Available options are:
#' \itemize{
#'   \item "RNA": this is the scRNA-seq data.
#'   \item "spatial": this is the spatial data.}
#' @param save.graphs Default is False , determing whether to save the graph(cell-cell similarity network).
#'   seurat graph is different from scanpy graph. Their relationship are set {"distances": "knn", "connectivities": "snn"} roughly.
#' @importFrom hdf5r H5File h5attr
#' @importFrom Hmisc capitalize
#'
seurat_to_h5 <- function(seurat=NULL, h5=NULL, assay.name = NULL, save.graphs = FALSE, save.scale = FALSE){
  data = h5$create_group('data')
  var = h5$create_group('var')
  if(assay.name == 'spatial'){
    seurat_spatial_to_h5(data=seurat, h5 = h5,gr_name = assay.name)
    assay.name <- capitalize(assay.name)
  }
  if(assay.name %in% names(slot(object = seurat, name = 'assays'))){
    #--- data
    slot_assay <- slot(object = seurat, name = 'assays')[[assay.name]]
    if(length(slot_assay@scale.data)>0){
      if(save.scale){
        n1 = 'rawX'
        sdata = 'X'
      }else{
        n1 = 'X'
        sdata = NULL
      }
    }else{
      n1 = 'X'
      sdata = NULL
    }
    df_to_h5(df = slot(object = slot_assay, name = 'meta.features'), h5 = var, gr_name = n1)
    if(all(slot_assay@data[,1] == slot_assay@counts[,1])){
      matrix_to_h5(mat = slot(object = slot_assay, name = 'counts'), h5 = data, gr_name = n1)
    }else{
      layer = h5$create_group('layer')
      matrix_to_h5(mat = slot(object = slot_assay, name = 'counts'), h5 = layer, gr_name= 'counts')
      matrix_to_h5(mat = slot(object = slot_assay, name = 'data'), h5 = data, gr_name = n1)
    }
    if(!is.null(sdata)){
      matrix_to_h5(mat = slot(object = slot_assay, name = 'scale.data'), h5 = data, gr_name = sdata)
      var1 = slot(object = slot_assay, name = 'meta.features')
      df_to_h5(df = var1[rownames(slot(object = slot_assay, name = 'scale.data')), ], h5 = var, gr_name = sdata)
    }
    #--- save the cell annotation
    df_to_h5(df = slot(object = seurat, name = 'meta.data'), h5 = h5, gr_name = 'obs')
    #--- save the dimension reduction
    if(length(seurat@reductions)>0){
      dimR <- h5$create_group('dimR')
      dim.r <- slot(object = seurat, name = 'reductions')
      for(d in names(dim.r)){
        D = toupper(d)
        dimR[[D]] <- t(slot(dim.r[[d]],'cell.embeddings'))
      }
    }
    #--- save the graphs
    if(save.graphs){
      if(length(seurat@graphs)>0){
        graph_df <- slot(object = seurat, 'graphs')
        if(length(grep(assay.name, names(graph_df))) == 2){
          graphs <- h5$create_group('graphs')
          gra_list <- list()
          gra_list[[paste0(assay.name, '_nn')]] <- 'knn'
          gra_list[[paste0(assay.name, '_snn')]] <- 'snn'
          for(g in names(graph_df)){
            matrix_to_h5(mat=graph_df[[g]], h5 = graphs, gr_name = gra_list[[g]])
          }
        }
      }
    }
    #--- save the metadata colors
    if('uns' %in% names(seurat@misc)){
      uns <- h5$create_group('uns')
      for(colr in grep('colors', seurat@misc[['uns']], value = TRUE)){
        uns[[colr]] <- seurat@misc[['uns']][[colr]]
      }
    }
    if('layer' %in% names(seurat@misc)){
      if(!('layer' %in% names(h5))){
        layer = h5$create_group('layer')
      }
      for(l in names(seurat@misc[['layer']])){
        matrix_to_h5(mat = seurat@misc[['layer']][[l]], h5 = layer, gr_name = l)
      }
    }
  }
  else{
    stop('Please enter the correct assay name')
  }
  return('Done!')
}


#' create the VisiumV1 class inherting the package Seurat 'SpatialImage' class
#' @import Seurat
#' @importFrom methods new
VisiumV1 <- setClass(
  Class = 'VisiumV1',
  contains = 'SpatialImage',
  slots = list(
    'image' = 'array',

    'scale.factors' = 'scalefactors',
    'coordinates' = 'data.frame',
    'spot.radius' = 'numeric'
  )
)


#' The spatial meassage to the h5 file(Seurat)
#'
#' @param data The seruat object
#' @param h5 The h5 file
#' @param assay.name 'assay.name' must be the 'spatial' for saving the spatial message.
#' @importFrom hdf5r H5File h5attr
#'
seurat_spatial_to_h5 <- function(data, h5, gr_name = 'spatial'){
  h5spa <- h5$create_group(gr_name)
  for(sample_id in names(data@images)){
    sid_h5 <-  h5spa$create_group(sample_id)
    #--- save the image for the lowres
    sid_image_h5 <- sid_h5$create_group('image')
    sid_image_h5[['lowres']] <-  aperm(slot(data@images[[sample_id]], 'image'))
    #--- save the scale factor
    sid_scalefactors_h5 <- sid_h5$create_group('scalefactors')
    sf <- slot(data@images[[sample_id]], 'scale.factors')
    v1 <- c('spot_diameter_fullres','fiducial_diameter_fullres','tissue_hires_scalef', 'tissue_lowres_scalef')
    for(k in names(sf)){
      sid_scalefactors_h5[[grep(k,v1, value = TRUE)]] <- sf[[k]]
    }
    #--- save the coor
    coor_df = slot(data@images[[sample_id]], 'coordinates')
    colnames(coor_df) <- c('in_tissue', 'array_row','array_col', 'image_2','image_1')
    df_to_h5(df = coor_df[,c('in_tissue', 'array_row','array_col','image_1','image_2')], h5 = sid_h5, gr_name = 'coor')
  }
}


#' The h5 group spatial to the spatial message
#'
#' @param h5spa The h5 group for the spatial
#' @importFrom methods new slot
#' @return The spatial message
#'
h5_to_spatial <- function(h5spa){
  spatial_list<- list()
  for(sid in names(h5spa)){
    spatial_sid_list <- list()
    sid_h5 <- h5spa[[sid]]
    for(me in names(sid_h5)){
      if('image' == me){
        spatial_sid_list[[me]] <- aperm(sid_h5[[me]][['lowres']][,,])
      }
      if('scalefactors' == me){
        sf_list <- list()
        v1 = c("spot","fiducial","hires","lowres")
        for(sf in v1){
          sf_list[[sf]] <- sid_h5[[me]][[grep(sf,names(sid_h5[[me]]), value = TRUE)]][]
        }
        sf_o <- Seurat::scalefactors(spot = sf_list$spot, fiducial = sf_list$fiducial, hires = sf_list$hiresk, lowres = sf_list$lowres)
        spatial_sid_list[[me]] <- sf_o
      }
      if('coor' == me){
        coor_df <- h5_to_df(sid_h5[['coor']])
        coor_df <- coor_df[,c('in_tissue','array_row','array_col','image_2','image_1')]
        colnames(coor_df)<- c('tissue','row','col','imagerow','imagecol')
        spatial_sid_list[[me]] <- coor_df
      }
    }
    unnormalized.radius <- spatial_sid_list$scalefactors$fiducial * spatial_sid_list$scalefactors$lowres
    spot.radius <- unnormalized.radius/max(dim(x = spatial_sid_list$image))
    spatial_list[[sid]] <- new(Class = "VisiumV1", image = spatial_sid_list$image, scale.factors = spatial_sid_list$scalefactors,
                               coordinates = spatial_sid_list$coor,spot.radius = spot.radius)
  }
  return(spatial_list)
}






