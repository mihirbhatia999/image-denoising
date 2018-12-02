%loading image and it's dimensions
X = imread('lenna.tiff');
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
k = 
phi = zeros(n,k);
S = zeros(h,w);
X_2 = zeros(h,w);
Omega = zeros(k,512);