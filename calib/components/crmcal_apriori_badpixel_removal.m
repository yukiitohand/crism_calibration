function [ DN14c,BP ] = crmcal_apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,varargin )
% [ DN14c,BP ] = crmcal_apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,varargin )
%  Apply a priori bad pixel removal
%  Input parameters:
%    DN14b   : 14bit DN image (L,S,B) processed until quadratic ghost
%    BPdata1 : CRISMdata obj, CDR BP data, using DF before the measurement
%    BPdata2 : CRISMdata obj, CDR BP data, using DF after the measurement
%              can be empty
%  Output parameters
%    DN14c    : processed 14bit DN data
%    BP       : bad pixels
%
%  2019/11/06: YUKI ITOH: empty BPdata2 is supported.
%  *Detail*

interpOpt = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'INTERPOPT'
                interpOpt = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

if isempty(BPdata1.img),BPdata1.readimg(); end
if ~isempty(BPdata2)
    if isempty(BPdata2.img),BPdata2.readimg(); end
end
if isempty(DMdata.img),DMdata.readimg(); end

[L,S,B] = size(DN14b);

if ~isempty(BPdata2)
    BP = or(BPdata1.img,BPdata2.img);
else
    BP = logical(BPdata1.img);
end

DN14c = DN14b;
% for l=1:L
%     d = DN14b(l,:,:);
%     d(BP) = nan;
%     DN14c(l,:,:) = d;
% end

% inMask = DMdata.img~=1;
outMask = DMdata.img==1;

switch interpOpt
    case 1
        % interpolation
        for b = 1:B
            imb = DN14c(:,:,b)';
            imbip = imb;
            bp = squeeze(BP(1,:,b));
            bp = and((bp==1),squeeze(outMask(1,:,b)));
            gp_bool = and((bp==0),squeeze(outMask(1,:,b)));
            gp = find(gp_bool)';
            xq = 1:S; xq = xq(bp)';
            if length(gp) > 0.1*S
                imbip(bp,:) = interp1(gp,imb(gp_bool,:),xq,'linear','extrap');
                % imbip(l,bp) = y1ip;
            end
            DN14c(:,:,b) = imbip';
        end
%         for b = 1:B
%             imb = DN14c(:,:,b);
%             imbip = imb;
%             bp = squeeze(BP(1,:,b));
%             bp = and((bp==1),squeeze(outMask(1,:,b)));
%             gp = and((bp==0),squeeze(outMask(1,:,b)));
%             gp = find(gp);
%             xq = 1:S;
%             for l=1:L
%                 if length(gp) > 0.1*S
%                     y1ip = interp1(gp,imb(l,gp),xq(bp),'linear','extrap');
%                     imbip(l,bp) = y1ip;
%                 end
%             end
%             DN14c(:,:,b) = imbip;
%         end
    case 2
end

% restore the DMmask
% inMask = DMdata.img~=1;
% for l=1:L
%     d = DN14c(l,:,:);
%     d_ori = DN14b(l,:,:);
%     d(inMask) = d_ori(inMask);
%     DN14c(l,:,:) = d;
% end


end
