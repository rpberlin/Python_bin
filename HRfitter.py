function [bestX, bestHRstruct] = HRfitter(power,HRactual)
HRmin=90;
HRmax=200;
Pmax=350;
a1 = .0084;
a2 = .0040;
a3 = .0041;
a4 = .0049;
a5 = 19.8002;
nSamps = size(HRactual,1);
X0 = [HRmin; HRmax; Pmax; a1; a2; a3; a4; a5];
Xnorm = [.1;.1;.1;.0001;.0001;.0001;.00049;.01];
%Xnorm = [1;1;1;1;1;1;1;1];
delta = .000000001;
nParams = 8;


[err0,Hhat0]=HRerror(power,HRactual,X0(1),X0(2),X0(3),X0(4),X0(5),X0(6),X0(7),X0(8));
bestErr = err0
bestHhat = Hhat0;
bestX = X0;

% plot(1:nSamps,HRactual,1:nSamps,Hhat0);
nSearches = 300;
errHist = err0;
for j=1:nSearches

errVec=zeros(nParams,1);
for i=1:nParams
Xg=X0;
Xg(i)=Xg(i)+delta.*Xnorm(i);
[err,blah]=HRerror(power,HRactual,Xg(1),Xg(2),Xg(3),Xg(4),Xg(5),Xg(6),Xg(7),Xg(8));
errVec(i)=err-bestErr;
end
deltaVec0 = -1.0*errVec;

nLineSearchSteps = 20;
lineErr=zeros(nLineSearchSteps,1);
lineSearchStepSize0=.001;
lineSearchStepSize=lineSearchStepSize0.*(1.5.^(1:nLineSearchSteps)');
for i=1:nLineSearchSteps
    Xls=X0+lineSearchStepSize(i)*deltaVec0;
    [lineErr(i),blah]=HRerror(power,HRactual,Xls(1),Xls(2),Xls(3),Xls(4),Xls(5),Xls(6),Xls(7),Xls(8));
    if lineErr(i) < bestErr
        bestErr = lineErr(i);
        bestHhat = blah;
        bestX = Xls;
        besti=i;
    end
end
besti
bestErr
errHist = [errHist;bestErr];
X0=bestX;

end

subplot(1,2,1)
plot(1:nSamps,HRactual,1:nSamps,Hhat0,1:nSamps,bestHhat);
subplot(1,2,2)
plot(errHist)

bestHRstruct.HRmin=bestX(1);
bestHRstruct.HRmax=bestX(2);
bestHRstruct.Pmax=bestX(3);
bestHRstruct.a1=bestX(4);
bestHRstruct.a2=bestX(5);
bestHRstruct.a3=bestX(6);
bestHRstruct.a4=bestX(7);
bestHRstruct.a5=bestX(8);



end


function [error,Hhat] = HRerror(power,HRactual,Hmin, Hmax,Pmax,a1,a2,a3,a4,a5)


nPts=size(power,1);

Hhat=zeros(nPts,1);
X2=zeros(nPts,1);
Hhat(1)=HRactual(1);

for i=2:nPts
[Hdot, X2dot] = dHdtParametric(Hhat(i-1), Hmin, Hmax,power(i-1), Pmax,X2(i-1),a1,a2,a3,a4,a5);
Hhat(i)=Hhat(i-1)+Hdot;
X2(i) = X2(i-1)+X2dot;
end

error =sum((HRactual-Hhat).^2);

end


function [Hdot, X2dot] = dHdtParametric(H, Hmin, Hmax,P, Pmax,X2,a1,a2,a3,a4,a5)
% a1 = .0113;
% a2 = .0072;
% a3 = .0041;
% a4 = .0049;
% a5 = 19.8002;

X1 = (H-Hmin)/Hmax;
u = P/Pmax;
phi = a4*X1/(1+exp(-1*(X1-a5)));

X1dot = -a1*X1+a2*X2+a2*u;
X2dot = -a3*X2+phi;
Hdot = X1dot*Hmax;


end