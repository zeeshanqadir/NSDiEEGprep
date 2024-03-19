%% Computes the R-squared value for two one-dimensional distributions
%   Taken by HH from legacy DH/KJM code, 04/2021
function [rsq, p] = mnl_rsq(x1, x2)

    x1=double(x1(:));
    x2=double(x2(:));

    sum1=sum(x1);
    sum2=sum(x2);
    n1=length(x1);
    n2=length(x2);
    sumsqu1=sum(x1.*x1);
    sumsqu2=sum(x2.*x2);

    G=((sum1+sum2)^2)/(n1+n2);

    % /var(X,Y)
    rsq = (sum1^2/n1 + sum2^2/n2 - G)/(sumsqu1 + sumsqu2 - G);
    if mean(x2) > mean(x1) % inverse sign if second vector larger than 1st vector
       rsq = -rsq;
    end
    
    %% Calculate p-value from F-test
    
    SStot = (length(x1) + length(x2))*var([x1; x2], 1); % total sum of squares (var(_, 1) normalized by N)
    SSres = length(x1)*var(x1, 1) + length(x2)*var(x2, 1); % residual sum of squares using respective means
    F = (SStot - SSres)/(SSres/(n1+n2-2)); % F statistic using 1 df numerator and N-2 df denominator
    
    p = 1 - fcdf(F, 1, n1+n2-2); % p-value

end