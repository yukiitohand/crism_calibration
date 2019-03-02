function [ DN14g ] = multiply_by_integrationTime( RT14g,hkt )
% [ DN14g ] = mupltiply_by_integrationTime( RT14g,hkt )
%   Multip.y RT14g with integration time [milliseconds]

[L,S,Bands] = size(RT14g);

integ_t = [hkt.data.EXPOSURE;];
integ_t = integ_t(:);

rate = [hkt.data.RATE];
rate = rate(:);

[t] = get_integrationTime(integ_t,rate);

DN14g = RT14g .* repmat(t,[1,S,Bands]);

end