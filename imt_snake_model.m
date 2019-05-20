function [y1 y2] = imt_snake_model(f,x,y1,y2,option)
%% parse parameters
x_step = 10;
alpha1 = 0.2; % LII smoothness
alpha2 = 0.2; % MAI smoothness
mu = 1.4; % thickness uniformity
max_iter = 500;

if isstruct(option)
    arg_set = fieldnames(option);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=option.',arg_set{i},';']);
        disp(['Use ',arg_set{i},' = ',num2str(option.(arg_set{i}))]);
    end
end

if alpha1 <= 0, alpha1 = 1e-9; end
if alpha2 <= 0, alpha2 = 1e-9; end
if mu <= 0, mu = 1e-9; end

%% resampling
oldx = x;
y1 = interp1(x,y1,[x(1):x_step:x(end)]');
y2 = interp1(x,y2,[x(1):x_step:x(end)]');
x = [x(1):x_step:x(end)]';
if 1
    [fx,fy] = gradient(f);
    fy = fy./max(max(abs(fy))); % boundary force is normalized to 1
    fy1 = fy; % force field for LII
    fy2 = fy; % force field for MAI
else
    % LII EM and MAI EM
    p1 = 0.5; p2 = 0.3; % 0 < pi < 1
    Ip = I1([2:end,end],:);
    Im = I1([1,1:end-1],:);    
    h = fspecial('gaussian',5,1);    
    Ip = imfilter(Ip,h,'replicate');
    Im = imfilter(Im,h,'replicate');
    f1 = p1*Ip-(1-p1)*Im;
    f2 = p2*Ip-(1-p2)*Im;
    % force field
    [fx,fy1] = gradient(f1);
    fy1 = fy1./max(max(abs(fy1))); % boundary force is normalized to 1
    [fx,fy2] = gradient(f2);
    fy2 = fy2./max(max(abs(fy2))); % boundary force is normalized to 1
end

% compensate defficient columns
ncols = 1+x_step*ceil((size(f,2)-1)/x_step); 
fy1 = cat(2,fy1,repmat(fy1(:,end),[1,ncols-size(f,2)])); % add new columns
fy2 = cat(2,fy2,repmat(fy2(:,end),[1,ncols-size(f,2)])); % add new columns

x = [x;x(end)+x_step];
% N = length(x);
% A = generate_snake_matrix(alpha,beta,N+4);
% invAI = inv(A(3:end-2,3:end-2) + gamma * diag(ones(1,N)));
y1 = [y1;y1(end)];
y2 = [y2;y2(end)];
% for boundary condition
V_s = [x(1:4),ones(4,1)];
V_e = [x(end-3:end),ones(4,1)];

%% generate matrix
lambda1 = mu+alpha1;
lambda2 = mu+alpha2;
l = length(y1);
% construct A1
A1 = -diag(ones(l-1,1)*lambda1,-1);
A1 = A1 - diag(ones(l-1,1)*lambda1,1);
A1 = A1 + diag(ones(l,1)*(2*lambda1+1),0);
% construct A2
A2 = -diag(ones(l-1,1)*lambda2,-1);
A2 = A2 - diag(ones(l-1,1)*lambda2,1);
A2 = A2 + diag(ones(l,1)*(2*lambda2+1),0);
% construct B
B = -diag(ones(l,1)*2*mu);
B = B + diag(ones(l-1,1)*mu,-1);
B = B + diag(ones(l-1,1)*mu,1);
% construct E1 E2 F1 F2
E1 = inv(A1-B*inv(A2)*B);
E2 = inv(A2-B*inv(A1)*B);
F1 = inv(A2*inv(B)*A1-B); % F1 = E1*B*inv(A2);
F2 = inv(A1*inv(B)*A2-B); % F2 = E2*B*inv(A1);

% E1(abs(E1)<1e-4) = 0; E2(abs(E2)<1e-4) = 0; 
% F1(abs(F1)<1e-4) = 0; F2(abs(F2)<1e-4) = 0;

%% evolution
iter = 1;
while iter <= max_iter
    % boundary force
    f1_b = interp2([1:size(fy1,2)],[1:size(fy1,1)],fy1,x,y1,'*linear',0);
    f2_b = interp2([1:size(fy2,2)],[1:size(fy2,1)],fy2,x,y2,'*linear',0);
    % linear extrapolation for boundary condition
    c1 = V_s \ mean([y1(1:4),y2(1:4)],2);    
    c2 = V_e \ mean([y1(end-3:end),y2(end-3:end)],2);
    yc_0 = polyval(c1,x(1)-x_step);
    yc_l1 = polyval(c2,x(end)+x_step);
    t_0 = y2(2)-y1(2);
    t_l1 = y2(end-1)-y1(end-1);    
    y1_0 = yc_0-t_0/2;
    y2_0 = yc_0+t_0/2;
    y1_l1 = yc_l1-t_l1/2;
    y2_l1 = yc_l1+t_l1/2; 
    % iteration
    % implicit scheme        
    Fy1 = f1_b; Fy2 = f2_b;
    % handles boundary condition
    Fy1(1) = Fy1(1)+lambda1*y1_0-mu*y2_0;
    Fy1(end) = Fy1(end)+lambda1*y1_l1-mu*y2_l1;
    Fy2(1) = Fy2(1)+lambda2*y2_0-mu*y1_0;
    Fy2(end) = Fy2(end)+lambda2*y2_l1-mu*y1_l1;
    C = y1 + Fy1;
    D = y2 + Fy2;
    
    y1 = E1*C-F1*D;    
    y2 = E2*D-F2*C;
    
    iter = iter + 1;
end
%% reduce redundent data
y1 = interp1(x,y1,oldx);
y2 = interp1(x,y2,oldx);

% delete(lh1); delete(lh2);
% lh1 = line(x,y1,'Color','g');
% lh2 = line(x,y2,'Color','g');

