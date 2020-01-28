function [dqdt] = LNSE(t,q,L)

global calls

% com[ute time derivative
dqdt    = -L*q;

calls   = calls+1;

if calls==10||mod(calls,100)==0
   disp(['    #' num2str(calls) ' t=' num2str(t)]) 
end

