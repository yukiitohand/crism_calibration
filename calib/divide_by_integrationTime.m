function [ RT14g ] = divide_by_integrationTime( DN14g,hkt )
% [ RT14g ] = divide_by_integrationTime( DN14g,hkt )
%   Divide DN14g with integration time [milliseconds]

[L,S,Bands] = size(DN14g);

integ_t = [hkt.data.EXPOSURE;];
integ_t = integ_t(:);

rate = [hkt.data.RATE];
rate = rate(:);

[t] = get_integrationTime(integ_t,rate);

RT14g = DN14g ./ repmat(t,[1,S,Bands]);

end
