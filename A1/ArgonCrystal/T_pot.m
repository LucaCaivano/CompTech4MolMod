function [TT] = T_pot(q, Sigma, Epsilon)
%T_pot
%   Input: matrix of q positions
%          matrix of Sigma
%          array of Epsilon
%   Output: TT potential energy

    TT = 0;
    for ii = 1:length(Sigma)
        for jj = (ii+1):length(Sigma)
            TT = TT + 4*Epsilon(ii,jj)*...
                (...
                (Sigma(ii,jj)/norm(q(ii,:)-q(jj,:)))^12 - ...
                (Sigma(ii,jj)/norm(q(ii,:)-q(jj,:)))^6 ...
                );
        end
    end
end