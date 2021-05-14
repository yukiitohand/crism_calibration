function [RT14j,K] = crmcal_subtract_highorderlight(RT14i,LLdata)
% [RT14j,K] = crmcal_subtract_highorderlight(RT14i,LLdata)
%  The ** step of the calibration: subtract high order light
%  Input Parameters
%    RT14i  : 14bit RT [counts/milliseconds] image (L,S,B) 
%             background subtracted
%    LLdata : CRISMdata obj, CDR LL data
%  Output parameters
%    RT14j    : processed data
%  * Detail *
%   section 2.12 in "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf" and CDR LL
%   data.

if isempty(LLdata.img), LLdata.readimg(); end
if isempty(LLdata.ROWNUM_TABLE), LLdata.read_ROWNUM_TABLE(); end



[L,S,Bands] = size(RT14i);

K = zeros([L,S,Bands]);
for n=1:6
    crds = LLdata.img(2*n-1,:,:);
    cffs = LLdata.img(2*n,:,:);
%     tmp = repmat(cffs,[L,1,1]) .* RT14i;
    valids = squeeze(abs(cffs)>0);
    for s=1:S
        crd = squeeze(crds(1,s,valids(s,:)));
        %crd_index = arrayfun(@(x) find(LLdata.ROWNUM_TABLE==x),crd)+1;
        % I don't understand why +1 is necessary. 2018.11.24
        % +1 is removed.
        crd_index = arrayfun(@(x) find(LLdata.ROWNUM_TABLE==x),crd);
        tmp = repmat(cffs(1,s,valids(s,:)),[L,1,1]) .* RT14i(:,s,crd_index);
        tmp_isnan = isnan(tmp);
        tmp(tmp_isnan) = 0;
        if ~isempty(crd)
            K(:,s,valids(s,:)) = K(:,s,valids(s,:)) + tmp;
        end
    end
end

RT14j = RT14i - K;

end
    
    