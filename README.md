# mat_to_10Xh5
Function for converting sce files from scGEATool to 10X Genomics  HDF5 format, while preserving metadata.

``` matlab
mat_to_10Xh5(sce,'output.h5');
```
## Load output in Seurat(R)
```R
library(seurat)
so = Read10X_h5('output.h5')
so = CreateSeuratObject(CreateAssayObject(so))
meta = read.csv('meta_output.csv', row.names = 1)
so = AddMetaData(so,meta)
```

## Load output in Scanpy(python)
```python
import scanpy as sc
import pandas as pd
adata = sc.read_10x_h5('output.h5')
adata.obs = pd.read_csv('meta_output.csv',index_col=0)
```

## References
James J Cai, scGEAToolbox: a Matlab toolbox for single-cell RNA sequencing data analysis, Bioinformatics, Volume 36, Issue 6, March 2020, Pages 1948–1949, https://doi.org/10.1093/bioinformatics/btz830
