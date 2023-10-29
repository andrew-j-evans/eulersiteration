clear;
clc;
format longG;

R0 = [1131.34; -2282.343; 6672.423];
V0 = [-5.64305; 4.30333; 2.42879];
mu = 398600;

R0Mag = norm(R0);
V0Mag = norm(V0);
[a,e,i,RA,w,nu,mu] = RV2COE(R0,V0);

dtv  = [0.1;1;10;60;300];
T = 1; %(days)

ERF = zeros(3,length(dtv));
EVF = zeros(3,length(dtv));
RFPQW = zeros(3,length(dtv));
SMEF = zeros(1,length(dtv));

%Iteration Preallocation
for j = 1:length(dtv)
    dt = dtv(j);
    n=T*60*60*24/dt;
    X = [R0(1);R0(2);R0(3);V0(1);V0(2);V0(3)]; %State Vector
    fXt = [X(4,1);X(5,1);X(6,1); -mu/R0Mag^3 * X(1,1); -mu/R0Mag^3 * X(2,1); -mu/R0Mag^3 * X(3,1)];

   
    R = zeros(1,n);
    V = zeros(1,n);
    Rpqw = zeros(3,n);
    SME = zeros(1,n);

    R(1) = R0Mag;
    V(1) = V0Mag;
    Rpqw(:,1) = [X(1,1); X(2,1); 0];
    SME(1) = (V0Mag^2/2)-(mu/R0Mag);
    Xn = X;

    for t=1:(n-1)
        Xnp1 = Xn + fXt*dt;
        Xn = Xnp1;
        
        R(t+1) = norm(Xn(1:3));
        V(t+1) = norm(Xn(4:6));
        
        Rmag = norm(R(t+1));

        Rpqw(:,t+1) = [Xn(1,1); Xn(2,1); 0];
    
        fXt = [Xn(4,1);Xn(5,1);Xn(6,1); -mu/Rmag^3 * Xn(1,1); -mu/Rmag^3 * Xn(2,1); -mu/Rmag^3 * Xn(3,1)];
        SME(t+1) = (V(t+1)^2/2)-(mu/R(t+1));

    end


ERF(:,j) = Xn(1:3);
EVF(:,j) = Xn(4:6);

APE(:, j) = cross(ERF(:,j), EVF(:,j));
RFPQW(:,j) = Rpqw(:,n);
SMEF(j) = SME(n);

end

%% Tables and Outputs
format shortG;
dt = table(reshape(dtv',1,[]),'VariableNames',"dt");
disp(dt)

finalREul = table(ERF,'VariableNames',"Euler's Final R Vectors (km)");
finalVEul = table(EVF,'VariableNames',"Euler's Final V Vectors (km/s)");
disp(finalREul);
disp(finalVEul);

orbitalElements = table(a, e, i, RA, w, nu, 'VariableNames',["a (km)","e","i (deg)","RA (deg)","w (deg)","nu (deg)"]);
disp(orbitalElements);

eulerSME = table(SMEF,'VariableNames',"Specific Mechanical Energy for Euler's Method");
disp(eulerSME);

APETable = table(APE,'VariableNames',"Specific Angular Momentum for Euler's Method");
disp(APETable);


%% Functions

function [a,e,i,RA,w,nu,mu] = RV2COE(R,V)
    mu = 398600;

    %Magnitude
    r = norm(R);
    v = norm(V);
    
    % Determination of Classical Orbital Elements
    % 1) Semi Major Axis
        epsillon = (v*v) / 2 - mu / r;
        a = - mu / (2 * epsillon);

    % 2) Eccentricity
        RVDot = dot(R,V);
        ebar = (1 / mu) * (R*(v*v-mu/r) - RVDot*V);
        e = norm(ebar);
    
    % 3) True Anomaly
        hbar = cross(R,V);
        h = norm(hbar);
        i = acosd(hbar(3)/h);
   
   % 4) Right Ascention of the Ascending Node
        kbar = [0,0,1];
        nbar = cross(kbar,hbar);
        n = norm(nbar);
        RA = acosd(nbar(1)/n);

        if nbar(2) < 0
            RA = 360 - RA;
        end

   % 5) Argument of Perigee
        nedot = dot(nbar,ebar);
        w = acosd(nedot/(n*e));
        if ebar(3) < 0
            w = 360 - w;
        end

   % 6) True Anomaly

        eRdot = dot(ebar, R);
        nu = acosd(eRdot/(e*r));   
        
        if RVDot < 0
            nu = 360 - nu;
        end


end

function [R,V] = COE2RV(a,e,i,RA,w,nu,mu)
    if(i == 0 || e == 0)
        disp("Warning!!! This is a special case! i or e = 0!")
    end
    P = a*(1-e*e);
    r = P / (1 + e*cosd(nu));

    snu = sind(nu);
    cnu = cosd(nu);

    Rpqw = [r*cnu; r*snu; 0];
    mup = sqrt(mu / P);
    Vpqw = [mup*(-snu); mup*(e+cnu); 0];
    
    cRA = cosd(RA);
    sRA = sind(RA);
    cw = cosd(w);
    sw = sind(w);
    ci = cosd(i);
    si = sind(i);
 
    A = [cRA*cw - sRA*sw*ci, -cRA*sw-sRA*cw*ci, sRA*si;
         sRA*cw+cRA*sw*ci, -sRA*sw+cRA*cw*ci, -cRA*si;
         sw*si, cw*si, ci
    ];

    R = A * Rpqw;
    V = A * Vpqw;

end
