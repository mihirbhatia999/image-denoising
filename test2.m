% R_i 
X = imread('C:\Users\Test\Desktop\lena.tiff');
X = rgb2gray(X);
X = imresize(X, 0.125);
X=double(X);
[h,w] = size(X);

R  = Ri(8,h,w);

patch_no=64;
Q=zeros(h,w);

for i=1:64
    X=reshape(X',[],1);
    P=R(:,:,i)'*R(:,:,i);
    Q=Q+transpose(reshape(P*X,w,h));
    
end

sum(permute(R,[2,1,3])*R,3)

function [R] = Ri(patch_size,h,w)
    %image_patch = X(xmin:xmin+(patch_size-1),ymin:ymin+(patch_size-1));
    %image_patch = reshape(image_patch,[],1);
no_of_patches = (h*w)/(patch_size*patch_size);
R = zeros(patch_size*patch_size, h*w, no_of_patches);
for k = 1:no_of_patches
    c = mod(k,(w/(patch_size))) ;
    if c==0
        c=patch_size;
    end
    r = floor(k/(w/(patch_size))) + 1;
    if mod(k,(w/(patch_size)))==0 
        r=k/(w/(patch_size));
    end
    xmin = (r-1)*patch_size + 1; 
    ymin = (c-1)*patch_size + 1 ;
    
    %giving 1 
      for b = 1:patch_size 
       for a = 1:patch_size
           R((b - 1)*patch_size + a,(xmin-1 + b - 1)*w + ymin + (a-1),k) = 1 ;
       end 
      end
  
end 
end