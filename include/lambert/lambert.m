function [A,P,V1,V2,conv] = lambert(X1,X2,TOF,mu,JJ,tol,kmax)
% This routine implements Battin's Solution to Lambert's Problem.
%      Inputs:
%      X1      Vector to Pt#1 
%      X2      Vector to Pt#2 
%      TOF     Time of Flight 
%      mu      Gravitational parameter for central body
%      JJ      Integer which sets initial condition for x (input)
%               (set JJ = 1 for an ellipse) 
%               (set JJ = 0 for an parabola or hyperbola)(Boolean)
%      tol     Tolerance to exit iterations
%      kmax    Maximum Iterations
%      Outputs:
%      V1      Velocity Vector at Pt#1 
%      V2      Velocity Vector at Pt#2
%      conv    Boolean to indicate if convergence happened

% Set convergence flag to true, will be changed if doesn't converge
conv = true;

%  (FIND TA (RAD) -
%  USE DOT PRODUCT TO OBTAIN BETWEEN 0 AND 180, THEN
%  CHECK TO SEE IF BETWEEN 180 AND 360 USING Z-COMPONENT
%  OF CROSS-PRODUCT OF RY AND RZ
ta = acos(dot(X1,X2)/(norm(X1)*norm(X2)));
if X1(1)*X2(2)-X1(2)*X2(1) < 0
    ta=2*pi-ta;
end

%(FIND CHORD)
c = sqrt(norm(X1)^2 + norm(X2)^2 - 2*norm(X1)*norm(X2)*cos(ta));

%(FIND SEMI-PERIMETER)
s = (norm(X1)+norm(X2)+c)/2;

%(FIND LAMBDA)
lambda = (sqrt(norm(X1)*norm(X2))*cos(ta/2))/s;

%(FIND W)
w = atan((norm(X2)/norm(X1))^0.25)-(pi/4);

%(FIND L)
l = ((tan(2*w))^2+(sin(ta/4))^2)/((tan(2*w))^2 + (sin(ta/4))^2 + cos(ta/2));

%(FIND M)
M = (8*mu*(TOF^2))/((s^3)*((1+lambda)^6));

x = 1;
y = 0; %Need t0 define to compile
if ~JJ
    x = 0.1E-13;
end

V1 = zeros(3,1);
V2 = zeros(3,1);
DX = 9e99;
% for J=1:100
%     if SKIP
%     else
k = 0;
while abs(DX) > tol
    n = x/(((sqrt(1+x))+1)^2);
    PHI = (8*(sqrt(1+x)+1))/(3+(1/(5+n+((9/7)*n)/(1+(((16/63)*n)/ ...
        (1+(((25/99)*n)/(1+(((36/143)*n)/(1+(((49/195)*n)/ ...
        (1+(((64/255)*n)/(1+(((81/323)*n)/(1+(100/399)*n))))))))))))))));
    H1 = ((l+x)^2*(1+3*x+PHI))/((1+2*x+l)*(4*x+PHI*(3+x)));
    H2 = (M*(x-l+PHI))/((1+2*x+l)*(4*x+PHI*(3+x)));
    B = (27*H2)/(4*((1+H1)^3));
    if B <= -1
        B = -1; 
    end
    U = -B/(2*(sqrt(1+B)+1));
    KU = (1/3)/(1-(((4/27)*U)/(1-(((8/27)*U)/(1-(((2/9)*U)/ ...
        (1-(((22/81)*U)/(1-(((208/891)*U)/(1-(((418/1755)*U)/ ...
        (1-(((598/2295)*U)/(1-(700/2907)*U)))))))))))))));
    YNEW = ((1+H1)/3)*(2+(sqrt(1+B)/(1-(2*U*(KU^2)))));
    XNEW = sqrt(((1-l)/2)^2+(M/YNEW^2))-((1+l)/2);
    %(COMPARE XNEW WITH X TO ENSURE CONVERSION)
    DX = abs(XNEW - x);
    x = XNEW;
    y = YNEW;
    k = k+1;
    if k > kmax
        conv = false;
        break
    end
end
A = M*s*(1+lambda)^2/(8*x*y^2);
P = (2*norm(X1)*norm(X2)*y^2*(1+x)^2*(sin(ta/2)^2))/(M*s*(1+lambda)^2);
F = 1-(norm(X2)/P)*(1-cos(ta));
G = norm(X1)*norm(X2)*sin(ta)/sqrt(mu*P);
FDOT = sqrt(mu/P)*tan(ta/2)*((1-cos(ta))/P-1/norm(X2)-1/norm(X1));
GDOT = 1-(norm(X1)/P)*(1-cos(ta));
for K=1:3
    V1(K)=(1/G)*(X2(K)-(F*X1(K)));
    V2(K)=FDOT*X1(K)+GDOT*V1(K);
end
end