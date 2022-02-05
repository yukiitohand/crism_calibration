function [ DN14b,sumGHOST_DN ] = crmcal_remove_quadrantGhost( DN14a,GHdata,hkt,varargin )
% [ DN14b,sumGHOST_DN ] = crmcal_remove_quadrantGhost( DN14a,GHdata,hkt,varargin )
%  The third step of the calibration: remove detector quadrant electronics
%  ghost
%  Input parameters:
%    DN14a  : 14bit DN image (L,S,B) bias removed
%    GHdata : CRISMdata obj, CDR GH data
%    hkt    : housekeeping table data (from TRR3)
%  Output parameters
%    DN14b    : processed DN14a data
%
%  Optional Parameters
%    'BINNING' : ID of binning mode {0,1,2,3}
%                (default) 0
%    'BINX'    : binning size (PIXEL_AVERAGING_WIDTH in LBL)
%                (default) 1 
%     *Note: The relationship between BINNING and BINX
%        BINNING   BINX
%              0      1
%              1      2
%              2      5
%              3     10
%    Please specify only one of 'BINNING', and 'BINX'
%
%  * Detail *
%  DN14b = DN14a - sum(H(GH))
%

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = crism_get_binx(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});   
        end
    end
end


if isempty(GHdata.tab), GHdata.readTAB(); end
[L,S,Bands] = size(DN14a);

rate = [hkt.data.RATE];
A = crism_frame_rateTABformatter(rate,GHdata.tab,'IR_GHOST_A','BINX',binx);
B = crism_frame_rateTABformatter(rate,GHdata.tab,'IR_GHOST_B','BINX',binx);
C = crism_frame_rateTABformatter(rate,GHdata.tab,'IR_GHOST_C','BINX',binx);
D = crism_frame_rateTABformatter(rate,GHdata.tab,'IR_GHOST_D','BINX',binx);
E = crism_frame_rateTABformatter(rate,GHdata.tab,'IR_GHOST_E','BINX',binx);

% GHOST_DN = zeros([L,S,Bands]);
% for b=1:Bands
%     DN14a_b = DN14a(:,:,b);
%     GHOST_DN_b = A + D.*(B-C)./2.*log(cosh(E./D)) ...
%         + (B+C)./2.*DN14a_b-(D.*(B-C)./2.*log(cosh((DN14a_b-E)./D)));
%     GHOST_DN(:,:,b) = GHOST_DN_b;
% end

GHOST_DN = A + D.*(B-C)./2.*log(cosh(E./D)) ...
        + (B+C)./2.*DN14a-(D.*(B-C)./2.*log(cosh((DN14a-E)./D)));

% sum GHOST_DN
x = crism_get_quadrantPxl(1,'BINX',binx);
sumGHOST_DN = GHOST_DN(:,x,:);
for i=2:4
    x = crism_get_quadrantPxl(i,'BINX',binx);
    sumGHOST_DN = sumGHOST_DN + GHOST_DN(:,x,:);
end
sumGHOST_DN = repmat(sumGHOST_DN,[1,4,1]);

% if binx>1
%     h = fspecial('average',[binx,1]); % convolution vector
%     sumGHOST_DN_new = zeros([L,S/binx,B]);
%     for i=1:Bands
%         ghost_DN_b = sumGHOST_DN(:,:,b);
%         ghost_DN_b_c = conv2(ghost_DN_b,h,'valid'); % take convolution and 
%         sumGHOST_DN_new(:,:,b) = ghost_DN_b_c(1:binx:end,:); % downsample
%     end
% end

DN14b = DN14a - sumGHOST_DN;

end