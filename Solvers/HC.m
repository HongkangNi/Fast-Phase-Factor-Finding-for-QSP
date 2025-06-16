function phi = HC(coef, parity, eta)
% Transpose if coef is a row vector
if ~isvector(coef)
    error('invalid input')
end
flag = 0;
if isrow(coef)
    coef = coef.';
    flag = -1;
end
%% Weiss algorithm
bc = 1i*coef/2;
if parity == 0     
    bc(1) = bc(1)*2;   
end
d = length(bc) - 1;
N = ceil(d/eta);
N = max(N,2*d+2);
thd = 1;
while thd > 1e-12
    ext_bc = [bc;zeros(N-2*d-1-parity,1);bc(end:-1:2-parity)];
    bz = ifft(ext_bc)*N;
    if parity == 1
        bz = bz.*exp(1i*pi/N*(0:N-1)');  
    end
    logsqrt_b = log(sqrt(1-abs(bz).^2));
    r = fft(logsqrt_b)/N;
    r(2:N/2)=0;
    r(N/2+1:N) = 2*r(N/2+1:N);
    G = ifft(r)*N;
    coef_a = fft(exp(G))/N;
    thd = norm(coef_a(floor(N/4):ceil(N*3/4)))^2/norm(coef_a)^2;
    valuec = bz.*exp(-G);
    if parity == 1
        valuec = valuec.*exp(-1i*pi/N*(0:N-1)');
    end
    CN = fft(valuec)/N;
    C = CN(1:d+1);
    N = N * 10;
end

%% Half Cholesky
B_col = imag(C(end:-1:1));
b = half_chol([[1;zeros(d,1)],B_col],B_col);
phi_new2 = atan(b);
phi = phi_new2(end:-1:1);

if flag == -1
    phi = phi.';
end