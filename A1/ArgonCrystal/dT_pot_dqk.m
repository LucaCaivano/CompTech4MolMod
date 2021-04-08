function [dTT_dqk] = dT_pot_dqk(q, Sigma, Epsilon, k)
%dT_pot_dqk
%   Input: matrix of q positions
%          matrix of Sigma
%          array of Epsilon
%          index k of q_k
%   Output: dTT_dqk derivative of potential energy wrt q_k (array)

    dTT_dqk = [0 0];
    for jj = 1:length(Sigma)
        if (jj ~= k)
            dTT_dqk = dTT_dqk + 4*Epsilon(k,jj)*...
                (...
                -12*Sigma(k,jj)^12*norm(q(k,:)-q(jj,:))^(-14)*(q(k,:)-q(jj,:)) ...
                + 6*Sigma(k,jj)^6* norm(q(k,:)-q(jj,:))^(-8) *(q(k,:)-q(jj,:)) ...
                );
        end
    end
end