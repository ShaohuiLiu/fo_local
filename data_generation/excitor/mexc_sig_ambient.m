function f = mexc_sig_ambient(t,k)
% Syntax: f = mexc_sig(t,k)
% 1:20 PM 15/08/97
% defines modulation signal for exciter control
global exc_sig n_exc
f=0; %dummy variable
if n_exc~=0
%  exc_sig(:,k)=zeros(n_exc,1);
%  exc_sig(1,k)=0.1;
%end
pmech0 = [2.5 5.45 6.5 6.32 5.052 7.0 5.6 5.4 8.0 5 10.0 13.5 17.957 17.85 10.0 20.0];
 if t<=0
     exc_sig(:,k) = zeros(n_exc,1);
 elseif t >= 1 % apply ambient input
%     exc_sig(:,k) = zeros(n_exc,1);
    noise = sqrt(.001*(0.005) .* pmech0(1:n_exc)) * randn([n_exc,1]);
%     exc_sig(:,k) = .001 .* randn([n_exc,1]);
    exc_sig(:,k) = noise;
%     fprintf('Exciter fault %d, total %d.  \n',loc,n_exc)
 else
     exc_sig(:,k) = zeros(n_exc,1);
 end
end
return