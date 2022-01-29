clear all;
sigma2_dbm= -90;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

ct=1;
snrdbm = [0: 5: 30   ]; Kt=20;
R0 = 1.5; Rs=6; %U0's rate, and Us's data rate
al=2; %path loss
c0=3*10^8; fc = 300*10^9; rl = 5*exp(-3);%molecular absorption
G = 10^(35.35/10);%mainlobe antenna gain
R_d = 10; theta=3.17; %radius and the angle of the wedge
d0 = 5; h0PL = (c0/4/pi/fc)^2/(d0^(al)+1)*exp(-rl*d0); %U0's distance and path loss
lam=0.01;%user density
phid = 0.1; % blockage 
vari = 0.00001;%bisection search control parameter
pe =  1 - (1-2*qfunc(theta/sqrt(2*2^2)))/(1-2*qfunc(180/sqrt(2*2^2)));%misalignment probaiblity
tau = 1.1;% parameter for new NOMA
xi = 10^(-4); %sidelobe/mainlobe, used when misalignment happens
eta1 = (c0/4/pi/fc)^2; 
ep0 = 2^R0-1; eps = 2^Rs-1;
 
for k = 1:length(snrdbm)
    snr = 10^((snrdbm(k)-30)/10)/sigma_square;
    sum1 = 0;      sum2 = 0;     sum3 = 0;    sum4=0;sum5=0;azzz=0;
    %find the power coefficients by using bisection search
    alpha_search = 0;
    P0target = tau*(pe*(1-exp(-ep0/snr/h0PL/G/xi))  + (1-pe)*(1-exp(-ep0/snr/h0PL/G)));
    lb = 0; ub = 1/(1+ep0);
    alphas_est = (lb+ub)/2;
    diff = 1;
    while diff>vari
        x_est = ep0/snr/(1-(1+ep0)*alphas_est);
        P0realized = pe*(1-exp(-x_est/h0PL/G/xi))  + (1-pe)*(1-exp(-x_est/h0PL/G));
        if P0realized<P0target %
            lb = alphas_est;
        else
            ub = alphas_est;
        end
        old_est = alphas_est;
        alphas_est =  (lb+ub)/2;
        diff = abs(alphas_est-old_est);
    end
    ccc(k) = alphas_est;    
     
    K=5;
     %actually the best      
     %new scheme 
     alphas2 = alphas_est; alpha02 = 1- alphas2;
     eta2 = phid^2/theta/gammainc(R_d*phid,2);
     epx = max(eps/snr/alphas2, ep0/snr/(alpha02-ep0*alphas2));
     stepx = 4000;
     rx = [0:R_d/stepx:R_d];
     sumk = 0;
     for ir = 1 : length(rx)
         r = rx(ir);                 
         sumk = sumk + (1 - exp(- epx *...
             (1+r^al)*exp(rl*r)/G/eta1) ) * exp(-phid*r)*r* R_d/stepx;
     end
     pus(k) = eta2^K*theta^K*(sumk)^K; 
     
     pu0(k) = pe*(1-exp(-ep0/snr/h0PL/G/xi))  + (1-pe)*(1-exp(-ep0/snr/h0PL/G));
     temp1= ep0/snr/(alpha02-ep0*alphas2);
     pu0_new(k) = pe*(1-exp(-temp1/h0PL/G/xi))  + (1-pe)*(1-exp(-temp1/h0PL/G));
     
     %approximated 
     alphas_est = (tau-1)/tau/(1+ep0);
     alphas2 = alphas_est; alpha02 = 1- alphas2;
     eta2 = phid^2/theta/gammainc(R_d*phid,2);
     epx = max(eps/snr/alphas2, ep0/snr/(alpha02-ep0*alphas2));
     stepx = 4000;
     rx = [0:R_d/stepx:R_d];
     sumk = 0;
     for ir = 1 : length(rx)
         r = rx(ir);                 
         sumk = sumk + (1 - exp(- epx *...
             (1+r^al)*exp(rl*r)/G/eta1) ) * exp(-phid*r)*r* R_d/stepx;
     end
     pusx(k) = eta2^K*theta^K*(sumk)^K; 
     temp1= ep0/snr/(alpha02-ep0*alphas2);
     pu0_newx(k) = pe*(1-exp(-temp1/h0PL/G/xi))  + (1-pe)*(1-exp(-temp1/h0PL/G));
end 
%semilogy(snrdbm,psim, snrdbm,pa)%,snrdbm,psimx, snrdbm,pax,snrdbm,paxap)% 
semilogy( snrdbm,pu0,snrdbm,pu0_new,    snrdbm,pus, snrdbm,pu0_newx,    snrdbm,pusx )% 
%semilogy(  snrdbm,pu0_new,    snrdbm,pax )% 

%plot the subfigure from Mathworks
% ax1 = gca; % Store handle to axes 1.
% ax2 = axes('Position',[.7 .55 .2 .2])
%box on;
% hold(ax2, 'on')
%semilogy(ax2, snrdbm(5:7),pu0(5:7),snrdbm(5:7),pu0_new(5:7),snrdbm(5:7),pu0_newx(5:7)  )
%$U_{k^*}$, NOMA with $\alpha_s^*$
%$U_0$, NOMA with $\alpha_s^*$
