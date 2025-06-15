function phi = Weiss_INFFT(coef, eta)

%% Weiss algorithm
bc = coef/2;
bc(1) = bc(1)*2;   
d = length(bc) - 1;
N = d/eta;
thd = 1;
% tic
while thd > 1e-10
    ext_bc = [bc,zeros(1,N-2*d-1),bc(end:-1:2)];
    bz = ifft(ext_bc)*N;
    bz = bz.';
    logsqrt_b = log(sqrt(1-abs(bz).^2));
    clear ext_bc bz;
    r = fft(logsqrt_b)/N;
    r(2:N/2)=0;
    r(N/2+1:N) = 2*r(N/2+1:N);
    G = ifft(r)*N;
    clear logsqrt_b r;
    AN = fft(conj(exp(G)))/N;
    thd = norm(AN(floor(N/4):ceil(N*3/4)))^2/norm(AN)^2;
    N = N * 3;
end
% we = toc

%% superfast layer stripping
ac = real(AN(1:d+1)); % coefficient of a^*(z)
F2 = INFFT(ac,bc(end:-1:1)');
phi = atan(F2(end:-1:1));
