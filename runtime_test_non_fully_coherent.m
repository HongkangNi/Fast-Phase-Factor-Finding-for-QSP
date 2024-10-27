opts.print = 0;
opts.useReal = 0;
parity = 0;

d_init = 100;

n = 10;

runtime_Newton_complex= nan(n,1);
runtime_CM_complex= nan(n,1);
runtime_FPI_new= nan(n,1);
runtime_NLFT_old= nan(n,1);
runtime_HC= nan(n,1);
for kk = 1:n
    d = d_init*2^(kk-1);
    
    % -----------------------------------------
    % random non-fully-coherent example
    rng(1);
    coef = rand(1,d)-0.5;
    bc = 1i*coef/2;
    bc(1) = bc(1)*2;
    d = length(bc) - 1;
    N = 50*d;
    ext_bc = [bc,zeros(1,N-2*d-1),bc(end:-1:2)];
    bz = ifft(ext_bc)*N;
    bz = bz.';
    coef = coef/max(abs(bz))*0.5;
    % -----------------------------------------
    
    if kk <= 4
        % Riemann-Hilbert-Weiss algorithm
        tic
        bc = 1i*coef/2;
        bc(1) = bc(1)*2;   
        d = length(bc) - 1;
        N = 50*d;
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



        % Newton's method
        opts.method = 'Newton';
        [~,out_Newton] = QSP_solver(coef',parity,opts);
        runtime_Newton_complex(kk) = out_Newton.time;

        % Contraction mapping (FPI)
        opts.method = 'FPI';
        [~,out_FPI] = QSP_solver(coef',parity,opts);
        runtime_CM_complex(kk) = out_FPI.time;
    
    end
    
    % FFPI
    [~,~,~, runtime] = QSP_FFPI(coef',parity,opts);
    runtime_FPI_new(kk) = runtime;
    
    % Half Cholesky algorithm
    tic
    phi = HC(coef, 0.5);
    runtime_HC(kk) = toc;

    
 
end

%% plot figures

figure;
loglog(d_init*2.^(0:n-1),runtime_Newton_complex,'-o','MarkerSize',10,'LineWidth',2)
hold on
loglog(d_init*2.^(0:n-1),runtime_CM_complex,'-s','MarkerSize',10,'LineWidth',2)
hold on
loglog(d_init*2.^(0:n-1),runtime_FPI_new,'-d','MarkerSize',10,'LineWidth',2)
hold on
loglog(d_init*2.^(0:n-1),runtime_NLFT_old,'-^','MarkerSize',10,'LineWidth',2)
hold on
loglog(d_init*2.^(0:n-1),runtime_HC,'-*','MarkerSize',10,'LineWidth',2)

set(gca,'FontSize',15);
hold off
legend("Newton","FPI",...
    "FFPI","RHW","HC", ...
    'Interpreter', 'latex','FontSize',15, 'Location', 'SouthEast')
xlabel('$d$', 'Interpreter', 'latex','FontSize',22)
ylabel('runtime(s)', 'Interpreter', 'latex','FontSize',22)
hold off
