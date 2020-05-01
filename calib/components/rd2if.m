function [IoF] = rd2if(RDn,SFdata,r)
% [IoF] = rd2if(RDn,SFdata,r)
%   convert radiance to I/F
%  Input Parameters
%   RDn: radiance cube [L x S x B]
%   SFdata: CRISMdata obj, CDR SF data
%   r : distance from the Sun (AU)
%  Output Parameters
%   IoF : I/F cube [L x S x B]

if isempty(SFdata.img), SFdata.readimg(); end

[L,S,B] = size(RDn);

IoF = pi .* RDn ./ SFdata.img .* (r.^2);

% IoF = nan(size(RDn));
% for l=1:L
%     RDn_l = RDn(l,:,:);
%     IoF_l = pi .* RDn_l ./ SFdata.img .* (r.^2);
%     IoF(l,:,:) = IoF_l;
% end

end