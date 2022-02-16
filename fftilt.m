%%%%%%%%%%Numerical inverse laplace transform 
%
%Inputs:
%   F: function handle for frequeuncy domain soln
%   tm: max time
%   M: # of time steps
%   alfa: attenuation parameter
%Ouputs:
%   ft: fourier transformed soln
function[ft,t] = fftilt(F,tm,M,alfa)
er = 1e-6; %error tolerance
Q=2*M+2; %Have to double interval to limit the error
t=tm/M*(0:M); %time steps 
%t = linspace(0,tm,M);
NT =(t(2)-t(1))*Q; %maximum of theoretic time interval 
omega=2*pi/(NT); %frequency step 
c=alfa - log(er)/NT;
s=c-1i*omega*(0:Q-1);
Fsc=feval(F,s);
%%%%%%%%%%%Two dim
dim1=size(Fsc,1); 
dim3=size(Fsc,3); 
ft=fft(Fsc(:,1:Q,:),[],2); 
ft=ft(:,1:M+1,:); 
ft=2*real(ft)-repmat(real(Fsc(:,1,:)),[1,M+1,1]);
ft=repmat(exp(c*t)/NT,[dim1,1,dim3]).*ft;
ft(:,1,:)=2*ft(:,1,:);

end%function