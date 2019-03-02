function [img_bias] = minipipeline_calibration_IR_BI_woBPinterp(CDRBIdata,varargin)
% [img_bias] = minipipeline_calibration_IR_BI_woBPinterp(CDRBIdata,varargin)
%  re-calculate CDR BI data without performing bad pixel interpolation

%  Optional Parameters
%      'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                         access an only show the path, 0: nothing
%                         (default) 0
%      'OUT_FILE'       : path to the output file
%                         (default) ''
%      'Force'          : binary, whether or not to force performing
%                         pds_downloader. (default) false

dwld = 0;
force = false;
outfile = '';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case {'DWLD','DOWNLOAD'}
                dwld = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'OUT_FILE'
                outfile = varargin{i+1};
        end
    end
end

if isempty(CDRBIdata.basenamesCDR), CDRBIdata.load_basenamesCDR(); end
if isempty(CDRBIdata.basenames_SOURCE_OBS)
    CDRBIdata.load_basenames_SOURCE_OBS(); 
end


% get source BIdata download if not exists.
EDRBIdataList = [];
for i=1:length(CDRBIdata.basenames_SOURCE_OBS.BI)
    basename_edrbi = CDRBIdata.basenames_SOURCE_OBS.BI{i};
    [dirfullpath_local_obsbi,~,~,~,~]...
        = get_dirpath_observation(basename_edrbi,'Download',dwld,...
        'Force',force,'OUT_file',outfile);
    bidata = CRISMdata(basename_edrbi,dirfullpath_local_obsbi);
    EDRBIdataList = [EDRBIdataList bidata];
end

% read HKP
for i=1:length(EDRBIdataList)
    EDRBIdataList(i).download(dwld);
end


% select ones with the same frame rate.
frame_rate = CDRBIdata.lbl.MRO_FRAME_RATE{1};
EDRBIdataList_s = [];
for i=1:length(EDRBIdataList)
    if EDRBIdataList(i).lbl.MRO_FRAME_RATE{1} == frame_rate
        EDRBIdataList_s = [EDRBIdataList_s EDRBIdataList(i)];
    end
end

% get exposuretime
expo_timeList = zeros(length(EDRBIdataList_s),1);
for i=1:length(EDRBIdataList_s)
    integ = EDRBIdataList_s(i).lbl.MRO_EXPOSURE_PARAMETER;
    [t] = get_integrationTime(integ,frame_rate,'Hz');
    expo_timeList(i) = t;
end

% covert raw 12bit biases to 14bits, then mean
DN14a = [];
PPdata = CDRBIdata(i).readCDR('PP');
for i=1:length(EDRBIdataList_s)
    EDRBIdataList_s.read_ROWNUM_TABLE();
    DN12 = EDRBIdataList_s.readimg();
    [ DN14 ] = DN12toDN14( DN12,PPdata,EDRBIdataList_s.rownum_table );
    DN14a_i = robust_v2('mean',DN14,1,'NOutliers',2);
    DN14a = cat(3,DN14a,DN14a_i);
end

% compute the term with a0I
[L,S,B] = size(DN14a);
BSdata = CDRBIdata.readCDR('BS');
if isempty(BSdata.tab),BSdata.readTAB(); end
binx = CDRBIdata.lbl.PIXEL_AVERAGING_WIDTH;
a0I = rateQuadrantTABformatter(rate,BSdata.tab,'A0','BINX',binx);
row_lambdaList = reshape(rownum_table,[1 1 B])+1; % +1 is already performed.
term_a0IList = [];
for i=1:length(EDRBIdataList_s)
    integ_t = expo_timeList(i);
    term_integ_t = (502/480)*(480-integ_t);
    term_a0I = repmat(a0I,[1,1,B]) ...
    .*heaviside(repmat(row_lambdaList,[1,S,1])-repmat(term_integ_t,[1,S,B]));
    term_a0IList = cat(3,term_a0IList,term_a0I);
end

% perform linear regression??
Y3 = DN14a+term_a0IList;
c0 = nan(1,S,B);
c1 = nan(1,S,B);
T = [ones(L,1) expo_timeList];
for b=1:B
    Y2 = Y3(:,:,B);
    chat = T\Y2;
    c0(1,:,b) = chat(1,:);
    c1(1,:,b) = chat(2,:);
end

img_bias = c0;

end
