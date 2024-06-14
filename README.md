# mat_to_10Xh5
Function for converting sce files from scGEATool to 10X Genomics  HDF5 format, while preserving metadata.

``` matlab
mat_to_10Xh5(sce,'output.h5');
```
## Load output in R
```R
library(seurat)
so = Read10X_h5('output.h5')
so = CreateSeuratObject(CreateAssayObject(so))
meta = read.csv('meta_output.csv', row.names = 1)
so = AddMetaData(so,meta)
```

## Load output to python
```python
import scanpy as sc
import pandas as pd
adata = sc.read_10x_h5('output.h5')
adata.obs = pd.read_csv('meta_output.csv',index_col=0)
print(adata)
```
