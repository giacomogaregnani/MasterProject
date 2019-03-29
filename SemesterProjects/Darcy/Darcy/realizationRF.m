function [a,Nva]= realizationRF(Lmax,LC,nu,sigma,Nsim)
% Computation of the RANDOM FIELD

% Matern Covariance function
f=@(x,y) 1/(gamma(nu)*0.5*2^nu)*(sqrt(2*nu)*sqrt(x.^2+y.^2)/LC).^nu.*besselk(nu,sqrt(2*nu)*sqrt(x.^2+y.^2)/LC);

Nh=2^(Lmax+1);
s=5;
L=1;
if nu>0.5
    L=ceil(max(6*LC,1)); %LC=.3; Nsim=1000;
end
M=Nh*L;

x=linspace(-L,L,2*s*M+1); y=x;

[X,Y]=meshgrid(x(1:end-1),y(1:end-1));

%f=exp(-(abs(X)+abs(Y))/LC);
ff=f(X,Y);  ff(s*M+1,s*M+1)=1; f=ff;

ck=(-1).^[0:2*s*M-1]'*(-1).^[0:2*s*M-1].*fft2(f)/(2*s*M)^2;
ck=real(fftshift(ck));
ck=[ck,conj(ck(:,1))];
ck=[ck;conj(ck(1,:))];
ck=ck((s-1)*M+1:(s+1)*M+1,(s-1)*M+1:(s+1)*M+1);
cn=ck;
ck=ck(1:end-1,1:end-1);
%ck1=[ck,conj(ck(:,1))];
%ck1=[ck1;conj(ck1(1,:))];


y1=randn(M+1,M+1,Nsim);
y2=randn(M+1,M+1,Nsim);
y3=randn(M+1,M+1,Nsim);
y4=randn(M+1,M+1,Nsim);

y00=y1(1,1,:);

y1plus0=y1(2:end,1,:);
y2plus0=y2(2:end,1,:);
y3plus0=y3(2:end,1,:);
y4plus0=y4(2:end,1,:);

y10plus=y1(1,2:end,:);
y20plus=y2(1,2:end,:);
y30plus=y3(1,2:end,:);
y40plus=y4(1,2:end,:);

y10minus=y1(1,end:-1:2,:);
y20minus=y2(1,end:-1:2,:);
y30minus=y3(1,end:-1:2,:);
y40minus=y4(1,end:-1:2,:);

y1minus0=y1(end:-1:2,1,:);
y2minus0=y2(end:-1:2,1,:);
y3minus0=y3(end:-1:2,1,:);
y4minus0=y4(end:-1:2,1,:);

y1plusplus=y1(2:end,2:end,:);
y2plusplus=y2(2:end,2:end,:);
y3plusplus=y3(2:end,2:end,:);
y4plusplus=y4(2:end,2:end,:);

y1plusminus=y1(2:end,end:-1:2,:);
y2plusminus=y2(2:end,end:-1:2,:);
y3plusminus=y3(2:end,end:-1:2,:);
y4plusminus=y4(2:end,end:-1:2,:);

y1minusplus=y1(end:-1:2,2:end,:);
y2minusplus=y2(end:-1:2,2:end,:);
y3minusplus=y3(end:-1:2,2:end,:);
y4minusplus=y4(end:-1:2,2:end,:);

y1minusminus=y1(end:-1:2,end:-1:2,:);
y2minusminus=y2(end:-1:2,end:-1:2,:);
y3minusminus=y3(end:-1:2,end:-1:2,:);
y4minusminus=y4(end:-1:2,end:-1:2,:);



z=[(y1minusminus-y2minusminus+1i*y3minusminus+1i*y4minusminus)/2, (y1minus0+1i*y4minus0)/sqrt(2),(y1minusplus+y2minusplus-1i*y3minusplus+1i*y4minusplus)/2;...
   (y10minus+1i*y30minus)/sqrt(2), y00, (y10plus-1i*y30plus)/sqrt(2);...
   (y1plusminus+y2plusminus+1i*y3plusminus-1i*y4plusminus)/2, (y1plus0-1i*y4plus0)/sqrt(2),(y1plusplus-y2plusplus-1i*y3plusplus-1i*y4plusplus)/2];

Nva=1+2*2*M+4*M^2;

a=zeros(2*M,2*M,Nsim); %a1=zeros(2*Nh+1,2*Nh+1,Nsim);
%tic
%for i=1:Nsim
%a(:,:,i)=ifft2(ifftshift(exp(1i*[0:2*Nh-1]*pi)'*exp(1i*[0:2*Nh-1]*pi).*sqrt(ck).*z(1:end-1,1:end-1,i))*(2*Nh)^2);
%end
%toc
Z=zeros(2*M,2*M,Nsim); 
exp1i=(-1).^[0:2*M-1]'*(-1).^[0:2*M-1].*sqrt(ck);
for i=1:Nsim
Z(:,:,i)=exp1i.*z(1:end-1,1:end-1,i);
end
a=ifft2(ifftshift(Z)*(2*M)^2);


a=real(a);
a=[a,a(:,1,:)];
a=[a;a(1,:,:)];
a=sigma*a(M+1:M+1+Nh,M+1:M+Nh+1,:);

a=reshape(a,(Nh+1)^2,Nsim,1);
a=a';

a=exp(a);




end