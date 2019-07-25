function [result,labelImage2] = bte2(original,labelImage1)

sizeof = size(labelImage1);
sizeof = sizeof(1)*sizeof(2);
original1 = rgb2gray(original);

s = regionprops(labelImage1,'Solidity');
j = [s.Solidity];


x = original;
x = rgb2gray(x);
y = mean2(x);
w = im2double(x);
z = var(w);
z = sum(z);
a = y/z;
if 1.7 <= a <=2.72
    t = 195;
elseif 2.4999 <= a <= 2.0999
    t = 210;
else
    t = 200;
end

mean  = regionprops(labelImage1,original1,'MeanIntensity');
k = struct2cell(mean);
k = cell2mat(k);
f = find(k<=t);
labelImage2 =  labelImage1;
%p = find(k<=100);
%f = union(f,p);
f1 = numel(f);
  for u = 1:f1
       %for r = 1:(numel(labelImage1))
            %if labelImage1==f(u)
       q = find(labelImage1==f(u));
       labelImage2(q)=0;
  end
 labelImage2 = imclearborder(labelImage2,8);
 show = drawregionboundaries(labelImage2,original,[255 0 0]);
 %figure
 %imshow(show),title('Image after thresholding')
 %qq =find(labelImage2~=0);
 
l = find(labelImage2 ~= 0);
result = zeros(size(labelImage2));
result(l) = 1;
result = imclearborder(result,8);
 show = drawregionboundaries(result,original,[255 0 0]);
  %figure
  if max(result(:)) == 0
  imshow(show),title( 'No Detected Tumor','Color','white','FontSize',16 )
  else
 imshow(show),title( 'Detected Tumor','Color','white','FontSize',16)
  end
end
