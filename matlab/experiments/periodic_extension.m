function [E] = periodic_extension(Lp,N)
% This is just used to help determine the analytic form of the pinv.

if N <= 2*(Lp-1)
   error('Overlap not implemented');
end

E = [zeros(Lp-1) zeros(Lp-1,N-2*(Lp-1)) eye(Lp-1);
     eye(N);
     eye(Lp-1) zeros(Lp-1,N-2*(Lp-1)) zeros(Lp-1)];

end
