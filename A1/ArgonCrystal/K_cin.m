function [KK] = K_cin(p, M)
%K_cin
%   Input: array of p
%          array of masses M
%   Output: KK kinetic energy

    KK = 0;
    for ii = 1:length(p)
        KK = KK + 0.5 * p(ii)^2 / M(ii);
    end
end