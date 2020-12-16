function [x,mu,sigma] = rZscores( x )
    [n, p] = size(x);
    mu = nan(1, p);
    sigma = nan(1, p);            
    for featureIndex=1:p            
        [tmcd,smcd] = unimcd(x(~isnan(x(:, featureIndex)), featureIndex), ceil(n*0.5));
        mu(featureIndex) = tmcd;
        sigma(featureIndex) = smcd;
%        assert(smcd>1e-10);
    end    
    mask = (sigma<1e-12) | (isnan(sigma));
    sigma(mask) = 1;
    mu(isnan(mu)) = 0;
    x = (x - repmat(mu, n, 1)) ./ repmat(sigma, n, 1);
end

