function phi = HC(coef, eta)

%% Weiss algorithm
bc = 1i*coef/2;
bc(1) = bc(1)*2;   
d = length(bc) - 1;
N = d/eta;
thd = 1;
while thd > 1e-12
    ext_bc = [bc,zeros(1,N-2*d-1),bc(end:-1:2)];
    bz = ifft(ext_bc)*N;
    bz = bz.';
    logsqrt_b = log(sqrt(1-abs(bz).^2));
    r = fft(logsqrt_b)/N;
    r(2:N/2)=0;
    r(N/2+1:N) = 2*r(N/2+1:N);
    G = ifft(r)*N;
    coef_a = fft(exp(G))/N;
    thd = norm(coef_a(floor(N/4):ceil(N*3/4)))^2/norm(coef_a)^2;
    valuec = bz.*exp(-G);
    CN = fft(valuec)/N;
    C = CN(1:d+1);
    N = N * 10;
end

%% Half Cholesky
B_col = imag(C(end:-1:1));
b = half_chol([[1;zeros(d,1)],B_col],B_col);
phi_new2 = atan(b);
phi = phi_new2(end:-1:1);