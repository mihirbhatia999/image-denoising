% X = ones(100,100);
% [h,w] = size(X);
% %adding noise to our image 
% mean_noise = 0;
% var_noise = 0.01;
% Y = imnoise(X,'gaussian',mean_noise, var_noise);
% S = zeros(h,w);
% lambda_2 = 1;
% 
% cvx_begin
%         variable S(h,w)
%         minimize (power(2,norm(Y-X-S,'fro'))+(lambda_2*norm(S,1)))
%         %minimize (power(2,norm(Y-X-S,'fro')))
%         %minimize(lambda_2*norm(S,1))
%  cvx_end
%     

X = imread('lenna.tiff');
X = rgb2gray(X);
T = laplace(X);
