%% loading image and it's dimensions, X_clean is the clean image to test
X_clean = imread('lenna.tiff');
X_clean = rgb2gray(X_clean);
X_clean = imresize(X_clean, 0.125);
[h,w] = size(X_clean);

%% adding noise to our image 
mean_noise = 0;
var_noise = 0.01;
Y = imnoise(X_clean,'gaussian',mean_noise, var_noise);
imshowpair(Y,X_clean, 'montage')
    %initializing clean image with noisy image at first 
Y = double(Y);
X = double(imnoise(X_clean,'gaussian',mean_noise, 0.01));
%X = double(X_clean);
%{
Y = Noisy Image 
Penalty parameters = lambda_1, lambda_2, mu 
Smoothing parameter = delta 
Stopping tolerance = e1, e2 
Clean Image = X 
%}

%% from table 1 for p =1 
lambda_1 = 10; 
lambda_2 = 0.1;
sigma = 3 ; %for additive noise 
%mu = sigma/30;
mu = 0.001;
%delta = 0.12;  %p=1
delta = 0;
e1 = 0.001;
e2 = 0.001;

%initializing
%k > n ---> it limits the rank of dictionary to n
n = h;
k = 81;

%dictionary 
phi = zeros(power(n,0.5),power(k,0.5));
for i = 1:power(n,0.5)
    for j = 1:power(k,0.5)
      phi(i,j) = cos((i-1)*(j-1)*(pi/11)) ;
    end
end
size(phi)
phi = kron(phi,phi);

patch_size = 8;
S = zeros(h,w,'double');
X_2 = zeros(h,w);
X_1 = zeros(h,w);   %change

no_of_patches = (h*w)/(power(patch_size,2));
Omega = zeros(k,no_of_patches);

L = (1/9)*[-1 -1 -1; -1 8 -1; -1 -1 -1];
t = 0;
 a = 1;
j = 1;
delta_1 = e1 + 1;
delta_2 = e2 + 1;
tmax = 3;  % to change   
smax = 20 ; % change
p=2;
T=20;       %own value
s = 0 ;

while (delta_2>=e2) && (s<=smax)
    S_s1 = argmin_S(Y,X,S,lambda_2);  %EQ 17
    S_s1
    X_s1 = argmin_X(Y,X,S_s1,lambda_1,L,delta,mu,X_2,p);
    X_s1
    delta_2 = min(power(norm(X_s1-X,'fro'),2),power(norm(S_s1-S,'fro'),2));
    delta_2
    X = X_s1;
    S = S_s1;
    s=s+1;
end 

%% functions  
%Function argmin
function [S] = argmin_S(Y, X, S, lambda_2)
    fprintf('argmin S called')
    [h,w] = size(X);
    cvx_begin
        variable S(h,w)
        minimize (square_pos(norm(Y-X-S,'fro'))+(lambda_2*norm(S,1)))
        %minimize (power(2,norm(Y-X-S,'fro')))
        %minimize(lambda_2*norm(S,1))
    cvx_end
end
    
function [X] = argmin_X(Y,X,S,lambda_1,L,delta,mu,X_2,p)
    fprintf('argmin X function called')
    [h,w] = size(X);
    %kernel = -1 * ones(3);
    %kernel(2,2) = 8; 
    %LX = conv2(X, kernel, 'same');
     
    cvx_begin
        variable X(h,w)
        %minimize (power(2,norm(Y-X-S,'fro'))+(lambda_1*trace(power((p/2),((LX.')*LX)+(delta*delta*eye(w)))))+(mu*power(2,norm(X-X_2,'fro'))))
        minimize (square_pos(norm(Y-X-S,'fro'))+(lambda_1*(square_pos(norm((laplace(X)),'fro'))+(w*delta*delta)))+(mu*square_pos(norm(X-X_2,'fro'))))
    cvx_end
end 