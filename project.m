%loading image and it's dimensions
X = imread('\\Client\C$\Karthik\Big Data\lena.tiff');
X = rgb2gray(X);
[h,w] = size(X);

%adding noise to our image 
mean_noise = 0;
var_noise = 0.01;
Y = imnoise(X,'gaussian',mean_noise, var_noise);
imshowpair(Y,X, 'montage')

%{
Y = Noisy Image 
Penalty parameters = lambda_1, lambda_2, mu 
Smoothing parameter = delta 
Stopping tolerance = e1, e2 
Clean Image = X 
%}

%from table 1 for p =1 
lambda_1 = 10; 
lambda_2 = 0.1;
sigma = 3 ;%for additive noise 
mu = sigma/30;
delta = 0.12;
e1 = 0.001;
e2 = 0.001;

%initializing
%k > n ---> it limits the rank of dictionary to n
k = 600
phi = zeros(n,k);
S = zeros(h,w);
X_2 = zeros(h,w);
X_1 = zeros(h,w);   %change
Omega = zeros(k,512);
L = (1/9)*[-1 -1 -1; -1 8 -1; -1 -1 -1];
t = 0;
a = 1;
j = 1;
delta_1 = e1 + 1;
delta_2 = e2 + 1;
tmax = 10;  % to change
smax = 10 ; % change
p=1;
patch_size = 8;

%while loop
while (delta_1>=e1) && (t<=tmax)
    s=0;
    while (delta_2>=e2) && (s<=smax)
        S_s1 = argmin_S(Y,X,S,lambda_2);  %EQ 17
        X_s1 = argmin_X(Y,X,S,lambda_1,L,delta,mu,X_2,p);
        cvx_end
        delta_2 = min(power(norm(X_s1-X,'fro'),2),power(norm(S_s1-S,'fro'),2))
        X = X_s1;
        S = S_s1;
        s=s+1;
    end
    X_1_t1 = X;
    X_2_t1 = argmin(mu,X,X_1_t1,R_i,phi,Omega);
    delta_1 = min(power(norm(X_1_t1-X_1,'fro'),2),power(norm(X_2_t1-X_2,'fro'),2))
    X_1 = X_1_t1;
    X_2 = X_2_t1;
    t = t+1;
end

% line 16
X = X_2;

% SParse coding stage
% k_0 = sparsity value
Omega = OMP(X,phi,k_0);

% FUnction argmin
function [S] = argmin_S(Y, X, S, lambda_2)
    [h,w] = size(X);
    cvx_begin
        variable S(h,w)
        minimize (power(norm(Y-X-S,'fro'),2)+(lambda_2*norm(S,1)))
    cvx_end
end
    
function [X] = argmin_X(Y,X,S,lambda_1,L,delta,mu,X_2,p)
    [h,w] = size(X);
    cvx_begin
        variable X(h,w)
        minimize (power(norm(Y-X-S,'fro'),2)+(lambda_1*trace(power(((laplace(X).')*laplace(X))+(delta*delta*eye(w)),(p/2))))+(mu*power(norm(X-X_2,'fro'),2)))
    cvx_end
end 
function [LX] = laplace(X)
    kernel = -1 * ones(3);
    kernel(2,2) = 8; 
    LX = conv2(double(X), kernel, 'same');
end


function [X] = argmin(mu,X,X_1_t1,R_i,phi,Omega)
    [h,w] = size(X);
    cvx_begin
        variable X(h,w)
        minimize ((1/mu)*power(norm(X-X_1_t1,'fro'),2)+(lambda_1*trace(power(((laplace(X).')*laplace(X))+(delta*delta*eye(w)),(p/2))))+(mu*power(norm(X-X_2,'fro'),2)))
    cvx_end
end

function [image_patch] = patch(i,size, X)
    r = mod(i,(w/(size))) + 1;
    c = floor(i/(w/(size)));
    xmin = (c-1)*size + 1 
    ymin = (r-1)*size + 1 
    image_patch = X(xmin:xmin+7,ymin:ymin+7);
    
   
end