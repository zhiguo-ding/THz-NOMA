clear all;
sigma2_dbm= -90;%+10*log10(BW)+Nf; %Thermal noise in dBm
sigma_square=10^((sigma2_dbm-30)/10);

ct=1000000;
snrdbm = [0: 5: 30   ]; Kt=20;
R0 = 1.5; Rs=3; %U0's rate, and Us's data rate
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
    alphas_est = (tau-1)/tau/(1+ep0);
    
    for i =1 : ct           
        u  = rand;
        if u<pe
            Gx = G*xi;
        else
            Gx = G*1;
        end
        h0=complex(randn(1,1)*sqrt(0.5),randn(1,1)*sqrt(0.5))...
            *sqrt((c0/4/pi/fc)^2/(d0^(al)+1)*exp(-rl*d0)*Gx);  
        h0 = abs(h0)^2;

        
        % ppp
        muk = theta*lam*phid^(-2)*(gammainc(R_d*phid,2));%pi*R_d^2*lam;
        K = 5;%%poissrnd(muk);
        pppind = 1;
        a_pp=0; angle_u=0; clear hk;
        cx1=zeros(K,1); cy1=zeros(K,1); dx=zeros(K,1);
         while pppind-1<K                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             pppind is one larger than the actual number
             cx1(pppind,1) = sign(randn(1,1))*rand(1,1)*R_d;
             cy1(pppind,1) = sign(randn(1,1))*rand(1,1)*R_d;
             dx(pppind,1) = (cx1(pppind))^2+(cy1(pppind))^2;
             if dx(pppind,1)<R_d^2
                 angle1 = rand*2-1;
                 if abs(angle1)<theta/2 % only record the users in the theta
                     p = exp(-phid*sqrt(dx(pppind,1)));
                     u  = rand;
                     if u<p %thinning
                         gk=complex(randn(1,1)*sqrt(0.5),randn(1,1)*sqrt(0.5)); 
                         dist = dx(pppind,1)^(1/2);
                         hpl = eta1/(dist^(al)+1)*exp(-rl*dist);
                         hk(pppind) =  sqrt(hpl*G) * gk;  
                         pppind = pppind+1;
                     end
                 end
                                  
             end
         end         
         
         hbest = max(abs(hk).^2);
    
         %simulation CR-NOMA
         if log2(1+snr*h0)<R0 
             alphas = 0;             
         else
             alphas = (snr*h0-ep0)/snr/(1+ep0)/h0;
         end
         alpha0 = 1- alphas;
         
         if log2(1+snr*hbest*alpha0/(snr*hbest*alphas+1))<R0 
             sum1 = sum1 +1;
         elseif log2(1+snr*hbest*alphas)<Rs
             sum1 = sum1 +1;
         end  

         %simulation new CR-NOMA
         alphas2 = alphas_est;%1/4;
         alpha02 = 1- alphas2;
         if log2(1+snr*hbest*alpha02/(snr*hbest*alphas2+1))<R0
             sum2 = sum2 +1;
         elseif log2(1+snr*hbest*alphas2)<Rs
             sum2 = sum2 +1;
         end
    end
    
    psim(k) =sum1/ct;
    psimx(k) =sum2/ct;
    
     %analysis     
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
     pax(k) = eta2^K*theta^K*(sumk)^K;
     sumk = 0;
     for ir = 1 : length(rx)
         r = rx(ir);                 
         sumk = sumk + (1+r^al)*exp(rl*r) * exp(-phid*r)*r* R_d/stepx;
     end     
     paxap(k) = eta2^K*theta^K*epx^K/G^K/eta1^K/ (rl-phid)^(2*K)...
         *( (-rl+phid)^(-al)*(gamma(al+2) - gamma(al+2)*gammainc(R_d*(-rl+phid),al+2,'upper')  )...
         + exp(R_d*(rl-phid))*(R_d*(rl-phid)-1) +1)^K;
     %paxap(k)=eta2^K*theta^K/G^K/eta1^K *epx^K* (sumk)^K;
     
     %CR-NOMA      
     eta3 = (1-tau)*(ep0+2^R0*eps)/tau/snr/eps - 1/tau/snr;
     eta4 = (ep0+eps*2^R0)/snr;
     
     %term1           
     sum21 = pe*(1- exp(-ep0/snr /h0PL/G/xi)) ...
             + (1-pe)*(1- exp(-ep0/snr /h0PL/G));

     %term2
     sum22 = 0;
     y21 = [ep0/snr:(eta4-ep0/snr)/stepx:eta4];
     rx = [0:R_d/stepx:R_d];
     for iy = 1 : length(y21)
         sumk = 0;
         y = y21(iy);
         fy = pe/h0PL/G/xi*exp(-y/h0PL/G/xi) ...
             + (1-pe)/h0PL/G*exp(-y/h0PL/G) ;

         for ir = 1 : length(rx)
             r = rx(ir);                 
             sumk = sumk + (1 - exp(- ...
                 eps*(1+ep0)*y/(snr*(y-ep0/snr)) *...
                 (1+r^al)*exp(rl*r)/G/eta1) ) * exp(-phid*r)*r* R_d/stepx;
         end
         sum22 = sum22 + sumk^K*fy * (eta4-ep0/snr)/stepx;
     end

       %term3
     sum23=0;
     stepy = 0.000000001;
     y21 = [eta4:stepy:stepy*stepx*10];
     rx = [0:R_d/stepx:R_d];
     for iy = 1 : length(y21)
         sumk = 0;
         y = y21(iy);
         fy = pe/h0PL/G/xi*exp(-y/h0PL/G/xi) ...
             + (1-pe)/h0PL/G*exp(-y/h0PL/G) ;

         for ir = 1 : length(rx)
             r = rx(ir);                 
             sumk = sumk + ...
                 (1 - exp(- y*(1+r^al)*exp(rl*r)/G/eta1) ) * exp(-phid*r)*r* R_d/stepx;
         end
         sum23 = sum23 + sumk^K*fy * stepy;
     end
       
     pa(k) = sum21 + eta2^K*theta^K*(sum22 +  sum23); 
end 
%semilogy(snrdbm,psim, snrdbm,pa)%,snrdbm,psimx, snrdbm,pax,snrdbm,paxap)% 
semilogy(snrdbm,psim, snrdbm,pa,snrdbm,psimx, snrdbm,pax,snrdbm,min(1,paxap))% 