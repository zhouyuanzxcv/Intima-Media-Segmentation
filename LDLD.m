function [fh,x,y1,y2] = LDLD(varargin)
%% linked dual line detection
% Author: Zhou Yuan, 2011-2-27
%
% The code implements the following paper:
%
% Zhou, Y.; Cheng, X.; Xu, X. & Song, E., Dynamic programming in parallel boundary detection with application to ultrasound intima-media segmentation 
% Medical image analysis, Elsevier, 2013, 17, 892-906

addpath('..\');
figure_on = 1;
use_snake = 1;

lmax = 60; lambda = 4;
w1 = 0; w2 = pi/9;
dis_min = 6; dis_max = 20;

theta_num = 45; dir_tol = 1; % delta_rho = M/rho_num
epsilon = 0.2; t1 = 0.5;

%% parse input parameters
option = [];
if length(varargin) == 0
    I = imread('phantom.bmp');
end
if length(varargin) > 0
    if ischar(varargin{1})
        disp(['Process ',varargin{1}]);
        I = imread(varargin{1});
    elseif isa(varargin{1},'uint8')
        I = varargin{1};
    end
end
if length(varargin) > 1
    option = varargin{2};
end
if isstruct(option)
    arg_set = fieldnames(option);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=option.',arg_set{i},';']);
        disp(['Use ',arg_set{i},' = ',num2str(option.(arg_set{i}))]);
    end
end

if figure_on
    close all;
end

%% algorithm implementation
% construct edge map using directional operator
f = scale_multi(I);
% f = compute_edge_map(I);
f1 = f(:,:);

[M,N] = size(f);
rho_num = M; 
seg_num = ceil(N/lmax);

start_t = tic;
h = fspecial('gaussian',5,1);
I1 = imfilter(double(I),h,'replicate');
[Ix,Iy] = gradient(I1);
grad_dir = atan2(Iy,Ix);        
[y1,y2] = dld_dp(f1,grad_dir,seg_num,epsilon,t1,lambda,rho_num,theta_num,dir_tol,w1,w2,dis_min,dis_max);
% [y1,y2] = dld_dp_mex_no_w2(f1,grad_dir,seg_num,epsilon,t1,lambda,rho_num,theta_num,dir_tol,w1,dis_min,dis_max);
x = (1:size(f,2))';
elapsed_t = toc(start_t);
disp(['Elapsed time is ',num2str(elapsed_t),' seconds']);

tic
if use_snake
    [y1 y2] = imt_snake_model(f,x,y1,y2,option);
end
toc

global total_comp_t;
total_comp_t = total_comp_t + elapsed_t;

%% display results
fh = [];
if figure_on    
    fh = figure;
    if ~isempty(varargin) && ischar(varargin{1})
        [~, name, ~, ~] = fileparts(varargin{1});
        set(fh, 'name', name);
    end
    subplot(2,1,1); imshow(I); hold on; plot(x,y1); plot(x,y2);
    subplot(2,1,2); imshow(f); hold on; plot(x,y1); plot(x,y2);
end

