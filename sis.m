% Demo Sequential Importance Sampling for linear gaussian model
clear all
close all
clc

T=100;
sw=2;
alpha=0.5;

% simulate data
x=zeros(1,T);
y=zeros(1,T);
x=randn(1);
for k=2:T
    x(1,k)=alpha*x(1,k-1)+randn(1);
end
y=x+sw*randn(1,T);

N=1000;

% prior
xs1=zeros(T,N);
lw1=zeros(T,N);
w1=zeros(T,N);
wnorm1=zeros(T,N);
ess1=zeros(1,T);

% optimal
xs2=zeros(T,N);
lw2=zeros(T,N);
w2=zeros(T,N);
wnorm2=zeros(T,N);
ess2=zeros(1,T);

% SIS using prior
for k=1:T
    if (k==1)
        xs1(k,:)=randn(1,N);
        lw1(k,:)=-0.5*log(2*pi*sw^2).*ones(1,N)-0.5*(y(1,k)*ones(1,N)-xs1(k,:)).^2./(2*sw^2);
    else
        xs1(k,:)=alpha.*xs1(k-1,:)+randn(1,N);
        lw1(k,:)=lw1(k-1,:)-0.5*log(2*pi*sw^2).*ones(1,N)-0.5*(y(1,k)*ones(1,N)-xs1(k,:)).^2./(2*sw^2);
    end
    lmax=max(lw1(k,:));
    w1(k,:)=exp(lw1(k,:)-lmax);  % correct only up to a multiplicative factor for unnormalized weights
    wnorm1(k,:)=w1(k,:)./sum(w1(k,:));
    ess1(1,k)=1/sum(wnorm1(k,:).^2);
end

% SIS using optimal
ss2=sw^2/(sw^2+1);
ss=sqrt(ss2);
for k=1:T
    if (k==1)
        xs2(k,:)=ss2*y(1,k)*ones(1,N)/sw^2+ss.*randn(1,N);
        lw2(k,:)=-0.5*log(2*pi*(sw^2+1)).*ones(1,N)-0.5*(y(1,k)*ones(1,N)).^2./(2*(sw^2+1));
    else
        xs2(k,:)=ss2.*(alpha.*xs2(k-1,:)+y(1,k)/sw^2)+ss.*randn(1,N);
        lw2(k,:)=lw2(k-1,:)-0.5*log(2*pi*(sw^2+1)).*ones(1,N)-0.5*(y(1,k)*ones(1,N)-alpha.*xs2(k-1,:)).^2./(2*sw^2);
    end
    lmax=max(lw2(k,:));
    w2(k,:)=exp(lw2(k,:)-lmax);  % correct only up to a multiplicative factor for unnormalized weights
    wnorm2(k,:)=w2(k,:)./sum(w2(k,:));
    ess2(1,k)=1/sum(wnorm2(k,:).^2);
end

std(lw1(T,:))
std(lw2(T,:))

figure(1)
plot(ess1,'b')
hold on
plot (ess2,'r')