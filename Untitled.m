%% This script gnerates the Figure 2 in the paper
%% Hongkang Ni and Lexing Ying, Fast phase factor finding for quantum signal processing, arXiv preprint
%% arXiv:2410.06409 (2024).
%% https://arxiv.org/abs/2410.06409

addpath("Solvers");
opts.print = 0;
opts.useReal = 0;
parity = 0;

tau_init = 100;

n = 2;

runtime_Newton_complex= nan(n,1);
runtime_CM_complex= nan(n,1);
runtime_FPI_new= nan(n,1);
runtime_NLFT_old= nan(n,1);
runtime_HC= nan(n,1);
for kk = 1:n
    tau = ceil(tau_init*1.7^kk);
    targ = @(x) 0.999*cos(tau.*x);
    d = ceil(1.4*tau+log(1e14));
    f = chebfun(targ,d);
    coef = chebcoeffs(f);
    coef = coef(parity+1:2:end)';
    
    % Half Cholesky algorithm
    tic
    phi = HC(coef, 1e-3);
    runtime_HC(kk) = toc;
    
    if kk <= 5
        % Riemann-Hilbert-Weiss algorithm
        tic
        bc = 1i*coef/2;
        bc(1) = bc(1)*2;   
        d = length(bc) - 1;
        N = 1000*d;
        ext_bc = [bc,zeros(1,N-2*d-1),bc(end:-1:2)];
        bz = ifft(ext_bc)*N;
        bz = bz.';
        logsqrt_b = log(sqrt(1-abs(bz).^2));
        r = fft(logsqrt_b)/N;
        r(2:N/2)=0;
        r(N/2+1:N) = 2*r(N/2+1:N);
        G = ifft(r)*N;
        valuec = bz.*exp(-G);
        CN = fft(valuec)/N;
        C = CN(1:d+1);

        phi_old = zeros(d+1,1);
        for k = 0:d
            H = hankel(C(k+1:end));
            id = eye(d+1-k);
            ab = [id, -H; -H, id]\[1;zeros(2*(d+1-k)-1,1)];
            a0 = ab(1);
            b0 = ab(d+2-k);
            phi_old(k+1) = atan(-1i*b0/a0);
        end
        runtime_NLFT_old(kk) = toc;





        coef = coef';

        % Newton's method
        opts.method = 'Newton';
        [~,out_Newton] = QSP_solver(coef,parity,opts);
        runtime_Newton_complex(kk) = out_Newton.time;
        
    end
   
 
end

%% plot figures


figure;
loglog(ceil(tau_init*1.7.^(1:n)),runtime_Newton_complex,'o-','MarkerSize',10,'LineWidth',2)
hold on
loglog(ceil(tau_init*1.7.^(1:n)),runtime_NLFT_old,'^-','MarkerSize',10,'LineWidth',2)
hold on
loglog(ceil(tau_init*1.7.^(1:n)),runtime_HC,'*-','MarkerSize',10,'LineWidth',2)

set(gca,'FontSize',15);
hold off
legend("Newton","RHW","HC", ...
    'Interpreter', 'latex','FontSize',15, 'Location', 'SouthEast')
xlabel('$\tau$', 'Interpreter', 'latex','FontSize',22)
ylabel('runtime(s)', 'Interpreter', 'latex','FontSize',22)
hold off