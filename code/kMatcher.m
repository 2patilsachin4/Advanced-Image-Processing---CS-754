function [kindicesPerFrame] = kMatcher(allpatches,i,j,f,k)
% allpatches: vectorized patches for all frames - 64 X (R/4-1) X (C/4-1) X
% nFrames
% f: frame no.      k: no. of matches
[~,h,w,nFrames] = size(allpatches);
L1 = sum(abs(allpatches - allpatches(:,i,j,f)));
L1 = reshape(L1,[h*w nFrames]);
[~,Linear_indices] = mink(L1,k);
[idxR,idxC] = ind2sub([h w],Linear_indices);
kindicesPerFrame = [idxR(:),idxC(:),repelem(1:nFrames,k)'];
end

