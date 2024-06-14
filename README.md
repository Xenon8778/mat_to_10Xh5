# mat_to_10Xh5
Function for converting sce files from scGEATool to 10X Genomics  HDF5 format, while preserving metadata.

``` matlab
mat_to_10Xh5(sce,'myfile.h5');
```
## Load output in R
```R
library(seurat)
so = Read10X_h5('test.h5')
so = CreateSeuratObject(CreateAssayObject(so))
meta = read.csv('meta_test.csv', row.names = 1)
so = AddMetaData(so,meta)
```
