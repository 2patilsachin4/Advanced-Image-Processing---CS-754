function [Filtered] = AdaptiveMedianRAMF(Corrupted,wSizeMax)
    
    Corrupted = double(Corrupted);
    Filtered = zeros(size(Corrupted));
    Flag_I = zeros(size(Corrupted));
    % Perform level I tests till window size is equal to wSizeMax
    k = 3; % Start with window size of 3 X 3
    while(k <= wSizeMax)
        % ordfilt2 (and medfilt2) for finding min, max (and median) parallely for all pixels;
        % symmetric - for padding with same values as border pixels, this
        % will help in retaining the median, min & max values


        MinMtx = ordfilt2(Corrupted,1,ones(k,k),'symmetric'); % 1st element in the ordered pixel values in a neighborhood is the min val 
        MaxMtx = ordfilt2(Corrupted,k^2,ones(k,k),'symmetric'); % last element in the ordered pixel values in a neighborhood is the min val
        MedMtx = medfilt2(Corrupted,[k k],'symmetric');
        
        % Level I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T1 = (MedMtx-MinMtx>0) & (MedMtx-MaxMtx<0) & (1-Flag_I);

        % Level II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T2 = (Corrupted-MinMtx>0) & (Corrupted-MaxMtx<0); 
        Filtered = Filtered + Corrupted.*(T2 & T1) + MedMtx.*(~T2 & T1); % level B passes-fill with I(x,y) else fill with I_med(x,y)
%         Corrupted
        % Update window size and pixels that have already been filtered
        k = k+2;    
        Flag_I = Flag_I + T1; % The pixels that have passed test 1 and need no further processing are added to Flags_I
        
        if all(Flag_I(:))
            break;
        end
        
    end
    Filtered(~Flag_I) = MedMtx(~Flag_I); %In case some pixels neither satisfy Test 1 nor Test 2 and k = wSizeMax
%     sum(Flag_I(:),1)

end


