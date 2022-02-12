function [ RT14g ] = crmcal_divide_by_integrationTime( DN14g,hkt )
% [ RT14g ] = crmcal_divide_by_integrationTime( DN14g,hkt )
%   Divide DN14g with integration time [milliseconds]

[L,S,Bands] = size(DN14g);

integ_t = [hkt.data.EXPOSURE;];
integ_t = integ_t(:);
if L==1
    integ_t = unique(integ_t);
    if length(integ_t)>1
        error('Something wrong');
    end
end

rate = [hkt.data.RATE];
rate = rate(:);
if L==1
    rate = unique(rate);
    if length(rate)>1
        error('Something wrong');
    end
end

[t] = crism_get_integrationTime(integ_t,rate);

% RT14g = DN14g ./ repmat(t,[1,S,Bands]);
RT14g = DN14g ./ t;

end
