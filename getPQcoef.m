function [lp, lq] = getPQcoef(phi)

n = length(phi);
phi = reshape(phi,[1,n]);

if mod(n,2) ~= 0
    error("invalid phi!");
elseif n == 0
    lp = [1;0];
    lq = [0,0];
else
    phiodd = phi(1:2:end);
    phieven = phi(2:2:end);
    N = 2^(ceil(log2(n))-1);
    lp = [1;0;0;0]*ones(1,N);
    lq = zeros(4,N);
    lp(:,1:n/2) = exp(1i*phiodd).*[1i*sin(phieven);0.5*cos(phieven);zeros(1,n/2);0.5*cos(phieven)];
    lq(:,1:n/2) = exp(1i*phiodd).*[zeros(1,n/2);0.5*cos(phieven);zeros(1,n/2);-0.5*cos(phieven)];
    
    while N > 1
        [m,N] = size(lp);        
        flp1 = fft([lp(1:m/2,1:2:end); zeros(m,N/2); lp(m/2+1:end,1:2:end)]);
        flp2 = fft([lp(1:m/2,2:2:end); zeros(m,N/2); lp(m/2+1:end,2:2:end)]);
        flq1 = fft([lq(1:m/2,1:2:end); zeros(m,N/2); lq(m/2+1:end,1:2:end)]);
        flq2 = fft([lq(1:m/2,2:2:end); zeros(m,N/2); lq(m/2+1:end,2:2:end)]);
        lp = ifft(flp1.*flp2 - flq1.*conj(flq2));
        lq = ifft(flp1.*flq2 + flq1.*conj(flp2));               
        N = N/2;
    end
    lp = lp.';
    lq = lq.';
end