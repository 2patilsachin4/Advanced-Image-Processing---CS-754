function [Q,stop] = FixedPoint(P, Good_pixels,stdDev,tau,MaxIter,eps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n1 = size(P,1);
n2 = size(P,2); % no. of selected patches(K*nFrames)
Q = zeros(size(P));
sigma_ = (mean(stdDev.^2)).^0.5; % average of the variances of all elements ∈ Ω
p = sum(sum(Good_pixels,1),2)/size(Good_pixels(:),1); % ratio of the number of pixels in Ω over the total number of pixels in the patch matrix
mu = (sqrt(n1) + sqrt(n2))* sqrt(p)*sigma_;

for i=1:MaxIter
    P_omega = Q-P;
    P_omega(~Good_pixels) = 0;
    R = Q - tau*P_omega;
    [U,Sig,V] = svd(R,'econ'); % Dτ (X) = UΣτV'
    Sig = diag(Sig);
    Q_ = U*diag(max(Sig-(tau*mu),0))*V'; 

    if (norm(Q_ - Q,'fro') <= eps) % Frobenius norm
        Q = Q_;
        
        break;
    end 
    Q = Q_;
end
stop = i;
end

