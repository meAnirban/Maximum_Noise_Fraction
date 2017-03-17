clear all
image = imread('mandril_color.tif');
[row,col,band_num] = size(image);
% convert image pixels to floating numbers
image = double(image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     noise covariance        %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% adding a row to each band
for a=1:band_num
    DN(:,:,a)=[image(:,:,a);image(row,:,a)];
end
   

for band=1:band_num
    for c=1:row
        N(c,:,band)=DN(c,:,band)-DN(c+1,:,band);
    end
end


for n=1:band_num
    % finding bandwise mean of noise matrix
    mean_N(1,n) = mean2(N(:,:,n));
    % bandwise normalization of noise matrix
    X1(:,:,n)=image(:,:,n)-mean_N(1,n)*ones(row,col);
    Y1(:,:,n)=image(:,:,n)-mean_N(1,n)*ones(row,col);
end

% for covariance matrix
sum_N=0;
for t=1:band_num
    for s=1:band_num
         for rows=1:row
            for cols=1:col
                z_N = X1(rows,cols,t)*Y1(rows,cols,s);
                sum_N= sum_N+z_N;
            end
         end
         cov_N = sum_N/((row*col)-1);
         cov_mat_N(t,s)=cov_N;
         sum_N=0;
    end
end

% eigen value(val_N) and eigen vector(vect_N) of noise
[vect_N,val_N] = eig(cov_mat_N);
% columnise eigen value of N
val_N = diag(val_N);
% sorting eigen value in descending order with their indices
[sort_val_N,index_N]=sort(val_N,'descend');


% sorting eigen vectors according to corresponding sorted eigen values
for o=1:length(sort_val_N)
    sort_vect_N(:,o) = vect_N(:,index_N(o));
end


% transformation with noise eigenvector
for r1=1:row
    for c1=1:col
        for b1=1:band_num
            % pixel value of image
            norml_img_N(b1,1)= image(r1,c1,b1);
        end
        % transformed value of each pixel with noise vector
        % transformed = feture noise vector transpose * original image
        pct_N = sort_vect_N.'*norml_img_N;
        for count1=1:band_num
            % principal component (noise)
            pct_img_N(r1,c1,count1)=pct_N(count1,1);
        end
    end
end

% now pct_img_N divided by corresponding bandwise noise standard deviation
% for noise whitening

for l=1:band_num
    std_N = std2(N(:,:,l));
    % transformed image data followed by pca
    F(:,:,l)= pct_img_N(:,:,l)/std_N;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                          %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%       PCA                %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% mean and normalizing image
for i=1:band_num
    % finding bandwise mean 
    mean(1,i) = mean2(F(:,:,i));
    % bandwise normalization
    X(:,:,i)=F(:,:,i)-mean(1,i)*ones(row,col);
    Y(:,:,i)=F(:,:,i)-mean(1,i)*ones(row,col);
end

% for covariance matrix
sum1=0;
for m=1:band_num
    for k=1:band_num
         for ro=1:row
            for co=1:col
                z = X(ro,co,m)*Y(ro,co,k);
                sum1= sum1+z;
            end
         end
         cov = sum1/((row*col)-1);
         cov_mat(m,k)=cov;
         sum1=0;
    end
end
    
                

% eigen value(val) and eigen vector(vect)
[vect,val] = eig(cov_mat);
% columnise eigen value
val = diag(val);
% sorting eigen value in descending order with their indices
[sort_val,index]=sort(val,'descend');

% sorting eigen vectors according to corresponding sorted eigen values
for j=1:length(sort_val)
    sort_vect(:,j) = vect(:,index(j));
end

% transformation
for r=1:row
    for c=1:col
        for b=1:band_num
            % pixel value of normalized image
            norml_img(b,1)= X(r,c,b);
        end
        % transformed value of each pixel
        % transformed = feture vector transpose * normalized image
        pct1 = sort_vect.'*norml_img;
        for count=1:band_num
            % principal component
            pct_img(r,c,count)=pct1(count,1);
        end
    end
end


figure, imshow(pct_img(:,:,1),[]);
figure, imshow(pct_img(:,:,2),[]);
figure, imshow(pct_img(:,:,3),[]);