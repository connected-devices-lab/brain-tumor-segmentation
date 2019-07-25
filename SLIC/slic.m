% SLIC Simple Linear Iterative Clustering SuperPixels
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.
%              k - Number of desired superpixels.
%       seRadius - Regions morphologically smaller than this are merged with
%                  adjacent regions.
%         colopt - String 'mean' or 'median' indicating how the cluster
%                  colour centre should be computed. Defaults to 'mean'
%             mw - Optional median filtering window size. 
% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%             Sp - Superpixel attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%                   stdL  - Standard deviation of L* 
%                   stda  - Standard deviation of a* 
%                   stdb  - Standard deviation of b* 
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each
%                           superpixel.
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.

function [l, Am, Sp, d] = done(im, k, m, seRadius, colopt, mw, nItr, eim, We)
    

k = 300;
m = 10;
seRadius = 1.5;
colopt = 'mean';
gh=im;
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
    
    % Convert image to L*a*b* colourspace.  This gives us a colourspace that is
    % nominally perceptually uniform. This allows us to use the euclidean
    % distance between colour coordinates to measure differences between
    % colours.
    im = RGB2Lab(im); 

    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero
    if mw
        if length(mw) == 1
            mw(2) = mw(1);  % Use same filtering for L and chrominance
        end
        for n = 1:3
            im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);
        end
    end
    
    % Nominal spacing between grid elements assuming hexagonal grid
    S = sqrt(rows*cols / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row
    nodeCols = round(cols/S - 0.5);
    % Given an integer number of nodes per row recompute S
    S = cols/(nodeCols + 0.5); 

    % Get number of rows of nodes allowing 0.5 row margin top and bottom
    nodeRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/nodeRows;

    % Recompute k
    k = nodeRows * nodeCols;
    
    % Allocate memory and initialise clusters, labels and distances.
    C = zeros(6,k);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is No of pixels
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.
    
    % Initialise clusters on a hexagonal grid
    kk = 1;
    r = vSpacing/2;
    
    for ri = 1:nodeRows
        % The following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
        if mod(ri,2), c = S/2; else, c = S;  end
        
        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
            C(1:5, kk) = [squeeze(im(rr,cc,:)); cc; rr];
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
    
    % Now perform the clustering.
    S = round(S);
    
    for n = 1:nItr
       for kk = 1:k  % for each cluster

           % Get subimage around cluster
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rows); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)
           
           % Compute distances D between C(:,kk) and subimage
           if USEDIST
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its
           % previous value update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end 
       % Update cluster centres with mean values
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
    end
    
    % Cleanup small orphaned regions on each region using
    % morphological opening on each labeled region.
%     [l, Am] = mcleanupregions(l, seRadius);
    Am = l;

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;
        nm = sum(mask(:));
        if centre == MEANCENTRE     
            Sp(n).L = sum(L(mask))/nm;
            Sp(n).a = sum(A(mask))/nm;
            Sp(n).b = sum(B(mask))/nm;
            
        elseif centre == MEDIANCENTRE
            Sp(n).L = median(L(mask));
            Sp(n).a = median(A(mask));
            Sp(n).b = median(B(mask));
        end
        
        Sp(n).r = sum(Y(mask))/nm;
        Sp(n).c = sum(X(mask))/nm;
        
        % Compute standard deviations of the colour components of each super
        % pixel.
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.
         
    end
    
    bte2(gh,l);
%-- distance 1
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre

% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.
function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    
    
% distance2
%
% Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%            eim - Edge strength sub-image corresponding to im
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.
function D = dist2(C, im, r1, c1, S, m, eim, We)

    % Squared spatial distance
    %    ds is a fixed image
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    % Combine colour and spatial distance measure
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    % for every pixel in the subimage call improfile to the cluster centre
    % and use the largest value as the 'edge distance'
    rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image
    cCentre = C(4)-c1;
    de = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            v = improfile(eim,[c cCentre], [r rCentre]);
            de(r,c) = max(v);
        end
    end

    % Combine edge distance with weight, We with total Distance.
    D = D + We * de;

    
   

     
    
    
  
