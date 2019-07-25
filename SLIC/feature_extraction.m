% FEATURE EXTRACTION

% labelImage = featextract(original,labelImage,adjim,n)
% original = the original rgb image 
% labelImage = segmented rgb image
% labelImage is the region adjacency graph matrix
% n is the centroids of each segment(region) 
% the function takes an rgb image,its segmented image,region adjacency graph matrix
% and centroids and returns merged regions of minimum dissimilarity 

% marks out the tumour section if there is any in the image
% W = cellfun(@(m) mean(m(1, :)), z);
function [labelImage1,weightadjim] = featextract(original,labelImage,adjim)
sizeof = size(labelImage);
sizeof = sizeof(1)*sizeof(2);

%Initialisation for calculation of consistency of cues
original = rgb2gray(original);
lamda1 = 0.5;
lamda = 0.5;
n0 = 0;
N0 = 5;
phi = 0;
alpha = 0.05;
beta = 0.005;
A = log10((1-beta)/alpha);
B = log10((beta)*(1-alpha));

% Finding the pixels corresponding to every region of the segmented image
stats = regionprops(labelImage,original,'PixelValues');
q = struct2cell(stats);

% Finding the average color of each region
mean  = regionprops(labelImage,original,'MeanIntensity');
k = struct2cell(mean);
k = cell2mat(k);
% Finding the area of each region
area = regionprops(labelImage,'Area','BoundingBox');
area = [area.Area];
m = size(adjim);
v = size(area);
w = zeros(m(1),1);

% Calculate the degree of disimilarity between regions
for i = 1:m(1)-1
    a = adjim(i,:);
    %s= n(a(1),:);
    %s1= n(a(2),:);
    %w(i) = ((area(a(1))*area(a(2)))/(area(a(1))+area(a(2))))*((s1(1)-s(1))^2 + ((s1(2)-s(2))^2));
    w(i) = ((area(a(1))*area(a(2)))/(area(a(1))+area(a(2))))*((k(a(1))-k(a(2)))^2);
    %w(i) = abs((k(a(1))-k(a(2))));
end
weightadjim = [adjim w];
z = cov(labelImage);
zz = inv(z);
h = size(z);
h = h(1);
count = 0;
labelImage1 = zeros(size(labelImage));
%labelImage1 = labelImage;
g = adjim(:,1);
g = union(g,adjim(:,1));
g1 = numel(g);
% Checks if adjacent regions to be merged are consistent
for i = 1:g1
    s = find(adjim(:,1)==g(i));
    v = min(weightadjim(min(s):max(s),3),[],1);
    a = find(weightadjim(:,3) == v);
    a = weightadjim(a,:);
      
  if (abs((k(a(1))-k(a(2))))<10)
      D=1;
  else
      D=0;
  end    
  % If they are consistent the regions are merged
 if D == 1
     for ii = 1 : sizeof
         if (labelImage(ii)==a(1) || labelImage(ii)==a(2))
             p = find(labelImage == a(1));
             p1 = find(labelImage == a(2));
         end
     end
          if labelImage1(p) == 0;
              count = count + 1;
             labelImage1(p) = count;
             labelImage1(p1) = count;
          elseif labelImage1(p) ~= 0;
             labelImage1(p1) = min(labelImage1(p));
             
          end
            
             
    
 end
end
end





