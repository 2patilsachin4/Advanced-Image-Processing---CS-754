clc;
clear all;

%% Reading some frames from input video
[Video, VideoInfo] = yuv4mpeg2mov('carphone_qcif.y4m');

nFrames = 50; h = VideoInfo.height; w = VideoInfo.width;
m = 0; sigma = 10; s = 40/100; k = 5; %Different noise level parameters
SetOfFrames = zeros([h,w,3,nFrames],'uint8');

for i = 1:nFrames
      SetOfFrames(:,:,:,i) =  Video(i).cdata;
end

SetOfFrames = double(SetOfFrames);

%% Adding noise to original frames

% Poisson Noise with zero mean, variance k*noisy
np = (poissrnd(k*SetOfFrames) - k*SetOfFrames);

% Gaussian Noise with zero mean, variance sigma^2
ng = sigma*randn(h, w,3, nFrames);
CorruptedSet = SetOfFrames + np + ng;

% Using uniform normal distribution for Impulsive noise
X = rand(h, w,3, nFrames);
CorruptedSet(X < s) = 0;
CorruptedSet(X > 1 - s) = 255;


%%  Computation of denoised frame
W = 11; FrameNo = 16;
[R_denoised,R_Adap] = RobustDenoise(uint8(CorruptedSet(:,:,1,:)),h,w,nFrames,FrameNo,W);
disp('Red frame filtering done')
[G_denoised,G_Adap] = RobustDenoise(uint8(CorruptedSet(:,:,2,:)),h,w,nFrames,FrameNo,W);
disp('Green frame filtering done')
[B_denoised,B_Adap] = RobustDenoise(uint8(CorruptedSet(:,:,3,:)),h,w,nFrames,FrameNo,W);
disp('Blue frame filtering done')

%% PSNR Calculation
% PSNR for corrupted frame
[R_PSNR,~] = psnr(uint8(CorruptedSet(:,:,1,FrameNo)),uint8(SetOfFrames(:,:,1,FrameNo)));
[G_PSNR,~] = psnr(uint8(CorruptedSet(:,:,2,FrameNo)),uint8(SetOfFrames(:,:,2,FrameNo)));
[B_PSNR,~] = psnr(uint8(CorruptedSet(:,:,3,FrameNo)),uint8(SetOfFrames(:,:,3,FrameNo)));
fprintf("\nThe PSNR of corrupted image is %.3f\n",(R_PSNR+G_PSNR+B_PSNR)/3);

% PSNR for Adaptive median filtered frame
[R_PSNR,~] = psnr(uint8(R_Adap),uint8(SetOfFrames(:,:,1,FrameNo)));
[G_PSNR,~] = psnr(uint8(G_Adap),uint8(SetOfFrames(:,:,1,FrameNo)));
[B_PSNR,~] = psnr(uint8(B_Adap),uint8(SetOfFrames(:,:,1,FrameNo)));
fprintf("The PSNR of Adaptive median filtered image is %.3f\n",(R_PSNR+G_PSNR+B_PSNR)/3);

% PSNR for final filtered frame
[R_PSNR,~] = psnr(uint8(R_denoised),uint8(SetOfFrames(:,:,1,FrameNo)));
[G_PSNR,~] = psnr(uint8(G_denoised),uint8(SetOfFrames(:,:,1,FrameNo)));
[B_PSNR,~] = psnr(uint8(R_denoised),uint8(SetOfFrames(:,:,1,FrameNo)));
fprintf("The PSNR of final output is %.3f",(R_PSNR+G_PSNR+B_PSNR)/3)

%% Display comparison
AdapMed = cat(3,R_Adap,G_Adap,B_Adap);
% clear R_Adap; clear G_Adap; clear B_Adap;
denoised_recon = cat(3,R_denoised,G_denoised,B_denoised);
% clear R_denoised; clear B_denoised; clear G_denoised;
output = cat(2,SetOfFrames(:,:,:,FrameNo),CorruptedSet(:,:,:,FrameNo),AdapMed,denoised_recon);
imshow(uint8(output))

