function f = scale_multi(I)
I1 = double(I); % image used for calculating force

s1 = 1; s2 = 15;

pw = 1:50; % possible widths
width1 = max(find(exp(-(pw.*pw)/(2*s1*s1))./(sqrt(2*pi)*s1)>0.0001));
width2 = max(find(exp(-(pw.*pw)/(2*s2*s2))./(sqrt(2*pi)*s2)>0.0001));

[x,y] = meshgrid(-width1:width1,-width1:width1);
dgau2D = -y.*exp(-(x.*x+y.*y)/(2*s1^2))/(2*pi*s1^2); % first derivative of gaussian
y = [-width2:width2]';
dgau1D = -y.*exp(-(y.*y)/(2*s2^2))/(sqrt(2*pi)*s2^3);     % the gaussian 1D filter

f1 = imfilter(I1, dgau2D, 'conv','replicate');
f1(f1<0) = 0;
f2 = imfilter(I1, dgau1D, 'conv','replicate');
f2(f2<0) = 0;
f = f1.*f2;

f = f./max(f(:));