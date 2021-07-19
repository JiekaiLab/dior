library(getopt)
library(dior)
spec <- matrix(
  c('readfile','r',1,'character','Reading the rds file',
    'targetobject', 't',1,'character','The single-cell data object which supprots Seurat and SingleCellExperiment'
  ),
  byrow = TRUE, ncol =5
)
opt <- getopt(spec)
# read the rds file
data <- readRDS(opt$readfile)
ind <- gregexpr('/', opt$readfile)
rdir <- substr(opt$readfile, start = 1, stop = max(ind[[1]]))
fname <- substr(opt$readfile, start = max(ind[[1]])+1, stop = nchar(opt$readfile))
wname <- gsub('.rds', '_tmp.h5', fname)
write_h5(data = data, object.type = opt$targetobject ,file = paste0(rdir, wname),
         assay.name = "RNA", save.graphs = TRUE, save.scale = FALSE)

