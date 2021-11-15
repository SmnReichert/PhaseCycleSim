% % % % % % % % % % % % % %
% % Helper class to calculate values like f_rs^q and dqqk
% % constants, etc
% % Methods are called by getValues.<Method>
 
classdef getValues < handle
properties(Constant)
    
    Na_gamma = 11.262*10^6; % Hz/T, gyromagnetic ratio for 23Na

end
methods(Static) 
    % get N random numbers in interval [a,b]
    function [r] = get_randomAB(a, b, N)
        r = a + (b-a) .* rand(N,1); 
    end
    
    % Calculate Lamor frequency for given B0 field strength
    function [w0] = get_w0(B0)
        w0 = getValues.Na_gamma * B0;
    end
       
    % Calculate p(tauC) for given tauCs and parameters for asymmetric log
    % Gauss distribution as specified in Rooney et al. (1991)
    function [ptauCs] = get_ptauC(tauCs, a, b, tcm)
        Norm = 2*a*b./(sqrt(pi)*(a+b)*tauCs);
        ptauCs =  Norm.*exp(-a^2*(log(tauCs/tcm)).^2) ;
        ptauCs(tauCs>tcm) = Norm(tauCs>tcm).*exp(-b^2*(log(tauCs(tauCs>tcm)/tcm)).^2) ;
    end    
    
    % Calculate wQ(tauC) for given tauCs
    % Step function with step at tauC0 and steps at wQ0 and wQ1
    function [wQs] = get_wQs(tauCs, tauC0, wQ0, wQ1)
        wQs = ones([ 1  length(tauCs)])*wQ0;
        wQs(tauCs>tauC0) = wQ1;
    end 
    
    % Calculate log Range of the Form 
    % a = startExp, b = endExp
    % 1*10^a, stepsize*2*10^a, ... stepsize*9*10^a, stepsize*1*10^(a+1), ... , stepsize*9*10^(b-1), 1*10^b
    function [logRange] = get_logRange(startExp, endExp, stepSize )
        r_tmp = (startExp:(endExp-1));
        t_tmp = 10.^r_tmp;
        logRange_tmp = [1:stepSize:9].'.*t_tmp;
        logRange = reshape(logRange_tmp,[1,length(1:stepSize:9)*length(r_tmp)]);
        logRange = [logRange 10^(endExp)];  
    end   
    
    % Calculate tauC and wQ for relaxation times T2s and T2f, dk20
    function [tauC, wQ] = get_tauC_wQ(T2f, T2s, w0)
        sizeF = length(T2f);
        sizeS = length(T2s);
        R11 = reshape(repelem(1./T2f, sizeS), [sizeS sizeF]);   % 1/s
        R12 = reshape(repelem(1./T2s, sizeF), [sizeF sizeS]).'; % 1/s
        a1 = R11./R12;
        b1 = R11-R12;
        b1(b1<0) = 0;
        a1(b1<0) = 0;
        b1(isnan(b1)) = 0;
        tauC = real(1./w0.*sqrt( 1./8 * ( -9 + 5.*a1 + sqrt( 25.*a1.^2 - 58.*a1 + 49 )))); % s
        x = (w0 * tauC).^2; % s
        wQ = sqrt(5*b1.* (1 + 4.*x)./(4*x*tauC)); % Hz
        wQ(isnan(wQ)) = 0;
    end
    
    % Calculate tauC and wQ for relaxation times T1s and T1f, dk20
    function [tauC, wQ] =  get_tauC_wQ_T1(T1f, T1s, w0)
        sizeF = length(T1f);
        sizeS = length(T1s);
        R01 = reshape(repelem(1./T1f, sizeS), [sizeS sizeF]);   % 1/s
        R02 = reshape(repelem(1./T1s, sizeF), [sizeF sizeS]).'; % 1/s
        a0 = R01./R02;
        b0 = R01-R02;
        b0(b0<0) = 0;
        a0(b0<0) = 0;
        % b1(isnan(b1)) = 0;
        tauC = real(1./w0.*sqrt( (a0-1)./(4-a0) )); % s
        x = (w0 * tauC).^2; % s
        wQ = sqrt(b0./(6/5*x.*tauC) .* (4.*x.^2+5*x + 1)); % Hz
        wQ(isnan(wQ)) = 0;
    end
    
    
    % Calculate relaxation Times T1s, T1f, T2s and T2f with given tauC and wQ 
    % based on dk20 definitions, 
    function [T1s, T1f, T2s, T2f] = get_relaxationTimes(tauC, wQ, w0)
        sizeTauC = length(tauC);
        sizeWQ = length(wQ);
        tauCs = reshape(repelem(tauC, sizeWQ), [sizeWQ sizeTauC]); % s
        wQs = reshape(repelem(wQ, sizeTauC), [sizeTauC sizeWQ]).'; % Hz
        
        x = (w0 .* tauCs).^2; % s
        a0 = (1+4*x)./(1+x);
        b0 = (6/5*x.*tauCs.*wQs.^2)./(4*x.^2+5*x+1);
        a1 = (2+9*x+4*x.^2)./(2+5*x);
        b1 = 1/5*tauCs.*wQs.^2 .* (4*x)./(1+4*x);
        
        T1f = (a0-1)./(a0.*b0);
        T1s = (a0-1)./(b0);
        T2f = (a1-1)./(a1.*b1);
        T2s = (a1-1)./(b1);
    end  
    
    % % get tauOpt for TQTPPI or IRTQTPPI, respectively 
    function [tauOpt] = get_tauOpt(Tif, Tis)
        tauOpt = log(Tis/Tif)/(1/Tif-1/(Tis));
    end 
    
    % spectral densities dk20
    % imaginary part Km can be absorbed by static Hamiltonian 
    % --> can be neglected in most situations
    function [Jm, Km] = get_J(m, tauC, wQ, w0)
        x = (w0 .* tauC).^2; % s
        Jm = (wQ.^2)/5.*tauC./(1+m.^2.*x); 
        Km = m.*w0.*tauC.*Jm;
    end 
    
    % Calculate relaxation rates for isotropic environment
    function [Rs] = get_Rs(q, tauC, wQ, w0)
        [J0, ~] =  getValues.get_J(0, tauC, wQ, w0);
        [J1, ~] =  getValues.get_J(1, tauC, wQ, w0);
        [J2, ~] =  getValues.get_J(2, tauC, wQ, w0);
        if q == 0
            R01 = 2*J1;
            R02 = 2*J2;
            R03 = 2*J1+2*J2;
            Rs = [R01 R02 R03];
        elseif q == 1 || q == -1
            R11 = J0+J1;
            R12 = J1+J2;
            R13 = J0+J1+2*J2;
            Rs = [R11 R12 R13];
        elseif q == 2 || q == -2
            R21 = J0+2*J1+J2;
            R22 = J0+J2;
            Rs = [R21 R22];
        elseif q == 3 || q == -3
            R31 = J1+J2;
            Rs = [R31];
        end
    end    
 
    % Calcluate relaxation rates for anisotropic environment
    function [Rs] = get_Rs_aniso(q, tauC, wQ, wQbar, w0)
        [J0, K0] =  getValues.get_J(0, tauC, wQ, w0);
        [J1, K1] =  getValues.get_J(1, tauC, wQ, w0);
        [J2, K2] =  getValues.get_J(2, tauC, wQ, w0);
        if q == 0
            R01 = 2*J1;
            R02 = 2*J2;
            R03 = 2*J1+2*J2;
            Rs = [R01 R02 R03];
        elseif q == 1 || q == -1
            R11 = J0+J1+J2-sqrt(J2^2-wQbar^2);
            R12 = J1+J2;
            R13 = J0+J1+J2+sqrt(J2^2-wQbar^2);
            Rs = [R11 R12 R13];
        elseif q == 2 || q == -2
            R21 = J0+J1+J2+sqrt(J1^2-wQbar^2);
            R22 = J0+J1+J2-sqrt(J1^2-wQbar^2);
            Rs = [R21 R22];
        elseif q == 3 || q == -3
            R31 = J1+J2;
            Rs = [R31];
        end
    end  
    
    
    % Calcluate Relaxation rates for anisotropic environmens and 
    % tauC distribution p(tauC)
    function [Rs] = get_Rs_aniso_ptauC(q, a, b, tcm, wQ, wQbar,  w0)
        tauCs = getValues.get_logRange(-14, -6, 0.1);
        ptauCs = getValues.get_ptauC(tauCs, a, b, tcm);
        deltaTauC = [tauCs(2)-tauCs(1) diff(tauCs)];
        wQs = wQ;
        % R = R(intgral J(tauC)dtauC)
        [J0s, K0] =  getValues.get_J(0, tauCs, wQs, w0);
        [J1s, K1] =  getValues.get_J(1, tauCs, wQs, w0);
        [J2s, K2] =  getValues.get_J(2, tauCs, wQs, w0);
        J0 =  sum(J0s.*ptauCs.*deltaTauC);
        J1 =  sum(J1s.*ptauCs.*deltaTauC);
        J2 =  sum(J2s.*ptauCs.*deltaTauC);
        if q == 0
            R01 = 2*J1;
            R02 = 2*J2;
            R03 = 2*J1+2*J2;
            Rs = [R01 R02 R03];
        elseif q == 1 || q == -1
            R11 = (J0+J1+J2-sqrt(J2^2-wQbar^2));
            R12 = (J1+J2);
            R13 = (J0+J1+J2+sqrt(J2^2-wQbar^2));
            Rs = [R11 R12 R13];
        elseif q == 2 || q == -2
            R21 = (J0+J1+J2+sqrt(J1^2-wQbar^2));
            R22 = (J0+J1+J2-sqrt(J1^2-wQbar^2));
            Rs = [R21 R22];
        elseif q == 3 || q == -3
            R31 = (J1+J2);
            Rs = [R31];
        end       
    end  
 

    % Calcluate Relaxation rates for anisotropic environment, 
    % tauC distribution p(tauC) and wQ distribution
    function [Rs] = get_Rs_aniso_ptauC_wQ(q, a, b, tcm, wQ0, wQ1, tauC0, wQbar,  w0)
        tauCs = getValues.get_logRange(-14, -6, 0.1);
        ptauCs = getValues.get_ptauC(tauCs, a, b, tcm);
        deltaTauC = [tauCs(2)-tauCs(1) diff(tauCs)];
        wQs = getValues.get_wQs(tauCs, tauC0, wQ0, wQ1);
        % R = R(intgral J(tauC)dtauC)
        [J0s, K0] =  getValues.get_J(0, tauCs, wQs, w0);
        [J1s, K1] =  getValues.get_J(1, tauCs, wQs, w0);
        [J2s, K2] =  getValues.get_J(2, tauCs, wQs, w0);
        J0 =  sum(J0s.*ptauCs.*deltaTauC);
        J1 =  sum(J1s.*ptauCs.*deltaTauC);
        J2 =  sum(J2s.*ptauCs.*deltaTauC);
        if q == 0
            R01 = 2*J1;
            R02 = 2*J2;
            R03 = 2*J1+2*J2;
            Rs = [R01 R02 R03];
        elseif q == 1 || q == -1
            R11 = (J0+J1+J2-sqrt(J2^2-wQbar^2));
            R12 = (J1+J2);
            R13 = (J0+J1+J2+sqrt(J2^2-wQbar^2));
            Rs = [R11 R12 R13];
        elseif q == 2 || q == -2
            R21 = (J0+J1+J2+sqrt(J1^2-wQbar^2));
            R22 = (J0+J1+J2-sqrt(J1^2-wQbar^2));
            Rs = [R21 R22];
        elseif q == 3 || q == -3
            R31 = (J1+J2);
            Rs = [R31];
        end        
    end 
   
    % Relaxation functions f(t) for isotropic environment
    function [fqkks] = get_f(q, t, tauC, wQ, w0)
        [J0, K0] =  getValues.get_J(0, tauC, wQ, w0);
        [J1, K1] =  getValues.get_J(1, tauC, wQ, w0);
        [J2, K2] =  getValues.get_J(2, tauC, wQ, w0);
        if q == 0
            R01 = 2*J1;
            R02 = 2*J2;
            R03 = 2*J1+2*J2;
            f011 = 1/5*(exp(-R01*t)+4*exp(-R02*t));
            f013 = 2/5*(exp(-R01*t)-exp(-R02*t));
            f033 = 1/5*(4*exp(-R01*t)+exp(-R02*t));
            f022 = exp(-R03*t);
            fqkks = [f011; f013; f033; f022];
        elseif q == 1 || q == -1
            R11 = J0+J1;
            R12 = J1+J2;
            R13 = J0+J1+2*J2;
            f111 = 1/5*(3*exp(-R11*t)+2*exp(-R12*t));
            f113 = sqrt(6)/5*(exp(-R11*t)-exp(-R12*t));
            f133 = 1/5*(2*exp(-R11*t)+3*exp(-R12*t));
            f122 = exp(-R13*t);
            fqkks = [f111; f113; f133; f122];
        elseif q == 2 || q == -2
            R21 = J0+2*J1+J2;
            R22 = J0+J2;
            f222 = exp(-R21*t);
            f233 = exp(-R22*t);
            fqkks = [f222; f233];
        elseif q == 3 || q == -3
            R31 = J1+J2;
            f333 = exp(-R31*t);
            fqkks = [f333];
        end
    end
    
    % Relaxation fuinctions f(t) for anisotropic environment
    function [fqkks] = get_f_aniso(q, t, tauC, wQ, wQbar, w0)
        [J0, K0] =  getValues.get_J(0, tauC, wQ, w0);
        [J1, K1] =  getValues.get_J(1, tauC, wQ, w0);
        [J2, K2] =  getValues.get_J(2, tauC, wQ, w0);
        Rs = getValues.get_Rs_aniso(q, tauC, wQ, wQbar, w0);
        if q == 0 
            R01 = Rs(1);
            R02 = Rs(2);
            R03 = Rs(3);
            f011 = 1/5*(exp(-R01*t)+4*exp(-R02*t));
            f013 = 2/5*(exp(-R01*t)-exp(-R02*t));
            f033 = 1/5*(4*exp(-R01*t)+exp(-R02*t));
            f022 = exp(-R03*t);
            fqkks = [f011; f013; f033; f022];
        elseif q == 1 || q == -1
            R11 = Rs(1);
            R12 = Rs(2);
            R13 = Rs(3);           
            mu = J2./(sqrt(J2.^2-wQbar.^2));
            v  = wQbar./(sqrt(J2.^2-wQbar.^2));
            f111 = 1/5*(3/2*(1+mu).*exp(-R11*t)+2*exp(-R12*t)+3/2*(1-mu).*exp(-R13*t));
            f122 = 1/2*((1-mu).*exp(-R11*t)+(1+mu).*exp(-R13*t));
            f133 = 1/5*((1+mu).*exp(-R11*t)+3*exp(-R12*t)+(1-mu).*exp(-R11*t));
            f113 = sqrt(6)/5*(1/2*(1+mu).*exp(-R11*t)-exp(-R12*t)+1/2*(1-mu).*exp(-R13*t));
            f112 = 1i/2*sqrt(3)/5*v*(sign(q)*exp(-R11*t)-sign(q)*exp(-R13*t));
            f123 = 1i/sqrt(10)*v*(sign(q)*exp(-R11*t)-sign(q)*exp(-R13*t));       
            fqkks = [f111; f113; f133; f112; f123; f122];
        elseif q == 2 || q == -2 
            R21 = Rs(1);
            R22 = Rs(2);            
            mu = J1./(sqrt(J1.^2-wQbar.^2));
            v  = wQbar./(sqrt(J1.^2-wQbar.^2));            
            f222 = 1/2*((1+mu)*exp(-R21*t)+(1-mu)*exp(-R22*t));
            f233 = 1/2*((1+mu)*exp(-R21*t)+(1-mu)*exp(-R22*t));
            f223 = -1i/2*v*(sign(q)*exp(-R21*t)-sign(q)*exp(-R22*t));
            fqkks = [f222; f233; f223];
        elseif q == 3 || q == -3
            R31 = Rs(1);
            f333 = exp(-R31*t);
            fqkks = [f333];
        end
    end
    
    % Wigner d-matrix element d_mbar,m^j(theta)
    function [djmm] = get_WignerD(j, mbar, m, theta)
        smax = max([(m-mbar), (j-mbar), (j+m)]); % find max s value, such that all factorials are non negative
        smax = max([smax 0]);
        djmm = 0;
        for s = 0:smax
            if mbar-m+s < 0 || j+m-s < 0 || j-mbar-s < 0
                continue;
            end
            djmm = djmm + ((-1)^(mbar-m+s)/(factorial(j+m-s)*factorial(s)*factorial(mbar-m+s)*factorial(j-mbar-s)) * cosd(theta/2)^(2*j+m-mbar-2*s)*sind(theta/2)^(mbar-m+2*s) );
        end 
        djmm = djmm * sqrt(factorial(j+mbar)*factorial(j-mbar)*factorial(j+m)*factorial(j-m));
    end
    
    % % get the nxn pulse superoperator for a specific order m  
    function [pulseOperator] = get_pulseOperator(flipAngle, phase, m)
        ns = -3:3;
        pulseOperator = zeros([length(ns) length(ns)]);  
        for nIndex = 1:length(ns)
            n = ns(nIndex);
            for nBarIndex = 1:length(ns)
                nbar = ns(nBarIndex);
                if abs(n)<=m && abs(nbar)<=m
                    pulseOperator(nIndex, nBarIndex) = exp(1i *(n-nbar)* deg2rad(phase))*getValues.get_WignerD(m, n, nbar, flipAngle);
                end
            end
        end
    end

end
end