function f = mexc_sig_impulse_all(t,k,loc)
% Syntax: f = mexc_sig(t,k)
% 1:20 PM 15/08/97
% defines modulation signal for exciter control
global exc_sig n_exc
f=0; %dummy variable
if n_exc~=0
%  exc_sig(:,k)=zeros(n_exc,1);
%  exc_sig(1,k)=0.1;
%end
 if t<=0
     exc_sig(:,k) = zeros(n_exc,1);
 elseif t >= 1 && t <= 1.003 % apply fault
%     exc_sig(:,k) = zeros(n_exc,1);
    exc_sig(:,k) = 0.05;
    fprintf('Exciter fault %d, total %d.  \n',loc,n_exc)
 else
     exc_sig(:,k) = zeros(n_exc,1);
 end
end
return