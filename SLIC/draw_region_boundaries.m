% DRAWREGIONBOUNDARIES: Draw boundaries of labeled regions in an image
% Usage: maskim = drawregionboundaries(l, im, col)
% Arguments:
%            l - Labeled image of regions.
%           im - Optional image to overlay the region boundaries on.
%          col - Optional colour specification. Defaults to black.  Note that
%                the colour components are specified as values 0-255.
%
% Returns: 
%       maskim - If no image has been supplied maskim is a binary mask
%                image indicating where the region boundaries are.
%                If an image has been supplied maskim is the image with the
%                region boundaries overlaid 

function maskim = drawregionboundaries(l, im, col)
    
    % Form the mask by applying a sobel edge detector to the labeled image,
    % thresholding and then thinning the result.
    h = [-1 1];  
    gx = filter2(h ,l);
    gy = filter2(h',l);
    maskim = (gx.^2 + gy.^2) > 0;
    maskim = bwmorph(maskim, 'thin', Inf);
    
    % Zero out any mask values that may have been set around the edge of the
    % image.
    maskim(1,:) = 0; maskim(end,:) = 0;
    maskim(:,1) = 0; maskim(:,end) = 0;
    
    % If an image has been supplied apply the mask to the image and return it 
    if exist('im', 'var') 
        if ~exist('col', 'var'), col = 0; end
        maskim = maskimage(im, maskim, col);
    end