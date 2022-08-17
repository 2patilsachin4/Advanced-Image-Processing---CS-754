function [denoised_recon,Adap_op] = RobustDenoise(SetOfFrames,h,w,nFrames,FrameNo,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


FiltGray = zeros([h,w,nFrames],'uint8');
Subset1 = zeros([h w nFrames]);
for i=1:nFrames
    corr = SetOfFrames(:,:,i);
    FiltGray(:,:,i) = AdaptiveMedianRAMF(corr,W);
    %Finding pixels detected by Adaptive median filter
    Subset1(:,:,i) = (FiltGray(:,:,i) == uint8(corr));
end
Adap_op = FiltGray(:,:,FrameNo);
FiltGray = double(FiltGray);

%######################## Seems correct till here #######################################################################

%% Creating Vectorized patches
allpatches = zeros([64 (h/4 - 1) (w/4 -1) nFrames]);
misboolarr = false(size(allpatches));
%%
for f=1:nFrames
    k1 = 0;
    for i=1:4:h - 7
        k2 = 0;
        for j=1:4:w - 7
            
            patch = FiltGray(i:i+7,j:j+7,f);
            mis_patch = Subset1(i:i+7,j:j+7,f);
            allpatches(:,i-3*k1,j-3*k2,f) = patch(:);
            misboolarr(:,i-3*k1,j-3*k2,f) = mis_patch(:);
            k2 = k2+1;
        end
        k1 = k1+1;
    end
end
clear k1; clear k2;
%% Finding similar patches
denoisedpatchArr = zeros(size(allpatches), 'double');
K = 5;
selPats = K*nFrames;
tau = 1.5;MaxIter = 30;eps = 1e-5;
for f = FrameNo:FrameNo %change this to nFrames later%%%%%%%%%%%%%%%
    % for each patch
    for i = 1:(h/4 -1)
        for j = 1:(w/4 -1)
            %find K*nFrames patches for each reference patch%%%

            % Find the indices corresponding to K*nFrames best patches per ref patch
            kindicesPerFrame = kMatcher(allpatches,i,j,f,K);
            KSelectedPatches = zeros(64,size(kindicesPerFrame,1));
            GoodPatchPixels = false(64,size(kindicesPerFrame,1));

            % find the actual patches using indices
            for s = 1:selPats
            KSelectedPatches(:,s) = allpatches(:,kindicesPerFrame(s,1),kindicesPerFrame(s,2),kindicesPerFrame(s,3)); % [64 patch x y frame]
            GoodPatchPixels(:,s) = misboolarr(:,kindicesPerFrame(s,1),kindicesPerFrame(s,2),kindicesPerFrame(s,3));
            end
            
            Mean = sum(KSelectedPatches.*GoodPatchPixels,2)./sum(GoodPatchPixels,2);
            stdDev = sum(((KSelectedPatches-Mean).*GoodPatchPixels).^2,2)./sum(GoodPatchPixels,2);
            stdDev = stdDev.^0.5;
            GoodPatchPixels = GoodPatchPixels & (KSelectedPatches <= Mean + 2*stdDev) & (KSelectedPatches >= Mean - 2*stdDev);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%5problem here

            [Denoised_patch,~] = FixedPoint(KSelectedPatches,GoodPatchPixels,stdDev,tau,MaxIter,eps);
            % for finding the patch no. of the reference patch from K*nFrames (250) selected patches
            temp = (kindicesPerFrame == [i,j,f]);
            patchNo = find(ismember(temp, [1 1 1],'rows'));
%              patchNo = 1+5*(f-1); % Since, L1 norm is the least for ref patch
            

            denoisedpatchArr(:,i,j,f) = Denoised_patch(:,patchNo);
        end 
    end
end

%%
%Averaged Reconstruction
denoised_recon = zeros([h,w]);
weighted_array= zeros(size(denoised_recon));
for f=FrameNo:FrameNo
    k1 = 0;
    for i = 1:4: h -7
        k2 = 0;
        for j=1:4:w - 7
            denoised_recon(i:i+7,j:j+7) = denoised_recon(i:i+7,j:j+7) + reshape(denoisedpatchArr(:,i-3*k1,j-3*k2,f), [8 8]);
            weighted_array(i:i+7,j:j+7) = weighted_array(i:i+7,j:j+7) + 1;
            k2 = k2+1;
        end
       k1 = k1+1;
    end
end

denoised_recon = denoised_recon./weighted_array;
end

