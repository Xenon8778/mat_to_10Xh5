function mat_to_10Xh5(sce, filename)
% MAT_TO_10XH5 Write and sce from workspace as a 10X HDF5 dataset. 
% H5WRITE(SCE,FILENAME) writes to an entire dataset and 
% a separate '.csv' file containing Cluster IDs, Cell Types and
% Experimental Batch IDs for each cell type
% Example:  Write to an entire dataset.
%       mat_to_10Xh5(sce,'myfile.h5');

if exist(filename, 'file')==2
  delete(filename);
end

if issparse(sce.X)
    disp('sparse matrix detected...');
else
    disp('converting to sparse...');
    sce.X = sparse(sce.X);
end

disp('converting to 10X sparse format...');
[indptr, indices, data] = convert_sparse_to_indptr(sce.X);

disp('writing matrix...');

h5create(filename, '/matrix/barcodes', size(sce.c_cell_id,1), Datatype='string');
h5write(filename, '/matrix/barcodes', sce.c_cell_id);

h5create(filename, '/matrix/data', size(data,1), Datatype='uint32');
h5write(filename, '/matrix/data', data);

h5create(filename, '/matrix/features/name', size(sce.g,1), Datatype='string');
h5write(filename, '/matrix/features/name', sce.g);

h5create(filename, '/matrix/features/id', size(sce.g,1), Datatype='string');
h5write(filename, '/matrix/features/id', sce.g);

feature_type = string(1:size(sce.g,1));
feature_type(:) = 'Gene Expression';
h5create(filename, '/matrix/features/feature_type', size(sce.g,1), Datatype='string');
h5write(filename, '/matrix/features/feature_type', feature_type);

genome = string(1:size(sce.g,1));
genome(:) = 'scRNA';
h5create(filename, '/matrix/features/genome', size(sce.g,1), Datatype='string');
h5write(filename, '/matrix/features/genome', genome);

h5create(filename, '/matrix/indices', size(indices,1), Datatype='uint32');
h5write(filename, '/matrix/indices', indices);

h5create(filename, '/matrix/indptr', size(indptr,2), Datatype='uint32');
h5write(filename, '/matrix/indptr', indptr);

h5create(filename, '/matrix/shape', size(size(sce.X),2), Datatype='uint64');
h5write(filename, '/matrix/shape', size(sce.X));

barcodes = sce.c_cell_id;
celltypes = sce.c_cell_type_tx;
clusters = sce.c_cluster_id;
batch = sce.c_batch_id;
meta_filename = append(filename(1:end-3),'_meta.csv');
metatable = table(barcodes,celltypes,clusters,batch);
writetable(metatable,meta_filename);
end

function [indptr, indices, data] = convert_sparse_to_indptr(X)

    % Check if X is sparse
    if ~issparse(X)
        error('Input matrix X must be a sparse matrix');
    end
    
    % Get matrix dimensions
    [~, n] = size(X);
    
    % Initialize indptr with 0 and n
    indptr = [0, n];
    
    % Find non-zero elements and their indices
    [row, col] = find(X);
     
    disp('sorting...')
    % Sort by columns for efficient construction
    [~, sort_idx] = sort(col);
    row = row(sort_idx);
    col = col(sort_idx);
    
    disp('accumulating...')
    % Accumulate column counts for indptr
    f = waitbar(0, 'Starting');
    for i = 1:n
        indptr(i+1) = indptr(i) + sum(col == i);
        waitbar(i/n, f, sprintf('Progress: %d %%', floor(i/n*100)));
    end
    close(f)
    
    % Assign indices and data. -1 to start indexing from 0
    indices = row-1;
    % data = full(X(row, col));  % Extract non-zero values
    data = nonzeros(X);

end