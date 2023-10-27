function f = mexc_sig_osci(t,k,mac_fault)
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
 elseif t >= 1 % apply osci input
%     exc_sig(:,k) = zeros(n_exc,1);
    noise = mac_fault(3) .* cos((2.*pi.*mac_fault(1)).*t + mac_fault(2));
    exc_sig(mac_fault(4),k) = noise;
%     disp('oscillation')
%     fprintf('Exciter fault %d, total %d.  \n',loc,n_exc)
 else
     exc_sig(:,k) = zeros(n_exc,1);
%      disp('0 input')
%      pause(0.5)
 end
end
return