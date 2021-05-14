function [RDn] = apply_DM(RDm,DMdata,varargin)
% [RSPl] = apply_DM(RDm,DM)
%  apply the detector mask defined in (CDR DM data).
%   Input
%     RDm   : radiance cube (L x S x B)
%     DMdata : CRISM obj, CDR DM data
%   Output Parameters
%     RDn: radiance cube on which non-scene pixels are replaced with
%          'MISSING_CONSTANT'
%   Optional Parameters
%     'MISSING_CONSTANT': value to be replaced at non-scene pixels
%                         (default) nan

missing_constant = nan;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MISSING_CONSTANT'
                missing_constant = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

[L,~,~] = size(RDm);

if isempty(DMdata.img), DMdata.readimg(); end;

mask = DMdata.img;

RDn = RDm;
for l=1:L
    RD_l = RDm(l,:,:);
    RD_l(mask~=1) = missing_constant;
    RDn(l,:,:) = RD_l;
end

end