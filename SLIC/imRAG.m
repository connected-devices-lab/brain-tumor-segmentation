function varargout = imRAG(img, varargin)
%IMRAG Region adjacency graph of a labeled image
%   Usage:
%   ADJ = imRAG(IMG);
%   computes region adjacencies graph of labeled 2D or 3D image IMG. 

%% Initialisations

% size of image
dim = size(img);

% number of dimensions
nd = length(dim);

% initialize array of neighbor regions
edges = [];

% Number of background pixels or voxels between two regions
% gap = 0 -> regions are contiguous
% gap = 1 -> there is a 1-pixel large line or surface between two adjacent
% 	pixels, for example the result of a watershed
gap = 1;
if ~isempty(varargin) && isnumeric(varargin{1})
    gap = varargin{1};
end
shift = gap + 1;

if nd == 2

    %% First direction of 2D image
    
    % identify transitions
    [i1 i2] = find(img(1:end-shift,:) ~= img((shift+1):end, :));
    
	% get values of consecutive changes
	val1 = img(sub2ind(dim, i1, i2));
	val2 = img(sub2ind(dim, i1+shift, i2));
	
    % keep only changes not involving background
    inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;
    edges = unique([val1(inds) val2(inds)], 'rows');
	

    %% Second direction of 2D image
    
    % identify transitions
    [i1 i2] = find(img(:, 1:end-shift) ~= img(:, (shift+1):end));
    
	% get values of consecutive changes
	val1 = img(sub2ind(dim, i1, i2));
	val2 = img(sub2ind(dim, i1, i2+shift));
    
    % keep only changes not involving background
    inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;
    edges = [edges; unique([val1(inds) val2(inds)], 'rows')];
    
    
elseif nd == 3
    %% First direction of 3D image
    
    % identify transitions
    [i1 i2 i3] = ind2sub(dim-[shift 0 0], ...
        find(img(1:end-shift,:,:) ~= img((shift+1):end,:,:)));
	
	% get values of consecutive changes
	val1 = img(sub2ind(dim, i1, i2, i3));
	val2 = img(sub2ind(dim, i1+shift, i2, i3));

    % keep only changes not involving background
    inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;
    edges = unique([val1(inds) val2(inds)], 'rows');
	

    %% Second direction of 3D image
    
    % identify transitions
    [i1 i2 i3] = ind2sub(dim-[0 shift 0], ...
        find(img(:,1:end-shift,:) ~= img(:,(shift+1):end,:)));
	
	% get values of consecutive changes
	val1 = img(sub2ind(dim, i1, i2, i3));
	val2 = img(sub2ind(dim, i1, i2+shift, i3));

    % keep only changes not involving background
    inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;
    edges = [edges; unique([val1(inds) val2(inds)], 'rows')];

    
    %% Third direction of 3D image
    
    % identify transitions
    [i1 i2 i3] = ind2sub(dim-[0 0 shift], ...
        find(img(:,:,1:end-shift) ~= img(:,:,(shift+1):end)));
	
	% get values of consecutive changes
	val1 = img(sub2ind(dim, i1, i2, i3));
    val2 = img(sub2ind(dim, i1, i2, i3+shift));
    
    % keep only changes not involving background
    inds = val1 ~= 0 & val2 ~= 0 & val1 ~= val2;
    edges = [edges; unique([val1(inds) val2(inds)], 'rows')];

end


% format output to have increasing order of n1,  n1<n2, and
% increasing order of n2 for n1=constant.
edges = sortrows(sort(edges, 2));

% remove eventual double edges
edges = unique(edges, 'rows');


%% Output processing

if nargout == 1
    varargout{1} = edges;
    
elseif nargout == 2
    % Also compute region centroids
    N = max(img(:));
    points = zeros(N, nd);
    labels = unique(img);
    labels(labels==0) = [];
    
    if nd == 2
        % compute 2D centroids
        for i = 1:length(labels)
            label = labels(i);
            [iy ix] = ind2sub(dim, find(img==label));
            points(label, 1) = mean(ix);
            points(label, 2) = mean(iy);
        end
    else
        % compute 3D centroids
        for i = 1:length(labels)
            label = labels(i);
            [iy ix iz] = ind2sub(dim, find(img==label));
            points(label, 1) = mean(ix);
            points(label, 2) = mean(iy);
            points(label, 3) = mean(iz);
        end
    end
    
    % setup output arguments
    varargout{1} = points;
    varargout{2} = edges;
end

