clc;
clear all;
%%

T = 5;
noise_std = 2;
noise_mean = 0;


%Read the video and extract first three frames
video = mmread('flame.avi');
frames = video.frames;
for i = 1:T
    gray_frames(:,:,i) = double(rgb2gray(frames(i).cdata));
end

[H,W,n_frames] = size(gray_frames);


%%
%Generate the coded single image using these
random_code = randi([0 1],H,W,T);
coded_snapshot = sum(gray_frames.*random_code,3);
noisy_snap_image = coded_snapshot + noise_std.*(randn(size(coded_snapshot)));


%%
%Generate the H*W*T coded vector

%For Cars Video cropping coordinates
%crop_frames = noisy_snap_image(179:274,127:342);    
%crop_random_code = random_code(179:274,127:342,:);

%For flames Video cropping ccoordinates
crop_frames = noisy_snap_image(130:257,56:287);
crop_random_code = random_code(130:257,56:287,:);

[new_H,new_W, n_frames] = size(crop_frames);

averaging_matrix = zeros(new_H,new_W,T);
reconstructed_img = zeros(new_H,new_W,T);

patch_size = 8;
si_matrix = kron(dctmtx(patch_size)',dctmtx(patch_size)');
%% 
%Main Algorithm to reconstruct the sub-frames using cropped noisy coded
%snapshot image
itr = 0;
for i = 1:new_H - 7
    for j = 1:new_W - 7
        A = [];
        patch = crop_frames(i:i+7,j:j+7)';
        coded_patch = crop_random_code(i:i+7,j:j+7,:);
        y = patch(:);
        
        for k = 1:T
            diag_vector_coded = coded_patch(:,:,k)';
            diag_patch_coded = diag(diag_vector_coded(:));
            A = [A diag_patch_coded*si_matrix];
        end
        
        theta = omp(A,y,noise_std^2);
        [tp2,~] = size(theta); 
        theta_3d = permute(reshape(theta',[1,tp2/T,T]),[2,1,3]);
        for m = 1:T
            patch_reconstructed_3d(:,:,m) = reshape(theta_3d(:,:,m),[patch_size,patch_size])';
            patch_reconstructed_3d(:,:,m) = idct2(patch_reconstructed_3d(:,:,m));
        end
        
        reconstructed_img(i:i+7,j:j+7,:) = reconstructed_img(i:i+7,j:j+7,:) + patch_reconstructed_3d;
        averaging_matrix(i:i+7,j:j+7,:) = averaging_matrix(i:i+7,j:j+7,:) + 1;
    end
end

%Average the reconstructed image over overlapping pixels
final_reconstructed_img = reconstructed_img./averaging_matrix;


%%
%RMSE calculation between original data and reconstructed data
%cropped_frames = gray_frames(179:274,127:342,:);    %For cars video
cropped_frames = gray_frames(130:257,56:287,:);     %For Flame video
rmse = sqrt(mean((final_reconstructed_img(:) - cropped_frames(:)).^2))
fprintf('The value of RMSE is: %f',rmse);

%%
%Display the original images alongside the reconstructed image
for i = 1:T
    figure;
    subplot(1,2,1), imshow(final_reconstructed_img(:,:,i)/255)
    title('Reconstructed img');
    subplot(1,2,2), imshow(cropped_frames(:,:,i)/255)
    title('Original img');
end


%%
%Store the coded snapshot image
%folder = '/pg/mtech/sachin/HW1-20220130T194005Z-001/HW1/results/';
%newimagename = [folder 'flame_coded_snapshot_5.jpg'];
%imwrite(uint8(255*(crop_frames/max(crop_frames(:)))),newimagename);


