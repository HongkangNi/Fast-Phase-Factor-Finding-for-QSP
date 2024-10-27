function cp = getUnitaryCoef(phi)

n = length(phi);
if mod(n,2) == 1
    [cp,~] = getPQcoef(phi(1:end-1));
    cp = cp*exp(1i*phi(end));
else
    [cp,cq] = getPQcoef(phi(1:end-2));
    cp = exp(1i*phi(end-1))*(cp + circshift(cp,-1))/2 + exp(-1i*phi(end-1))*(cq - circshift(cq,-1))/2;
    cp = cp*exp(1i*phi(end));
end
