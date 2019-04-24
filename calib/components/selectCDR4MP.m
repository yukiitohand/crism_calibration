function [SPdata,SSdata,SHdata] = selectCDR4MP(SPdata_input,varargin)

propSP_input = getProp_basenameCDR4(SPdata_input.basename);

% if the data is not VNIR, then raise error
if ~strcmpi(propSP_input.sensor_id,'S')
    error('The input must be VNIR data');
end
    

if propSP_input.binning==0 && propSP_input.wavelength_filter==0
    % if no binning and wavelength filter is applied, then just use the
    % input SPdata
    SPdata = SPdata_input; propSP = propSP_input;
elseif propSP_input.binning>0 && propSP_input.wavelength_filter==0
    % if binning is perfomed and wavelength filter is 0, then select SP
    % data with binning==0. If there is no such data, then use SPdata_input
    % directly.
    propSP_new = setfield(propSP_input,'binning',0);
    [dirpathSP_new,~,~,basenameSP_new] = get_dirpath_cdr_fromProp(propSP_new,varargin{:});
%     [basenamePtrnSP] = get_basenameCDR4_fromProp(propSP_new);
%     [dirpath_cdr,remote_subdir,basenameSP] = get_dirpath_cdr(SPdata_input.basename,varargin{:});
%     [basenameSP] = readDownloadBasename_v3(basenamePtrnSP,...
%                                     dirpath_cdr,remote_subdir,varargin{:});
    if ~isempty(basenameSP_new)
        SPdata = CRISMdata(basenameSP_new,dirpathSP_new);
        propSP = getProp_basenameCDR4(SPdata.basename);
    else
        SPdata = SPdata_input;
        propSP = propSP_input;
    end
else
    error('not implemented for wavelength_filter(%d)>0',propSP_input.wavelength_filter);
end

%VNIR SS data (as of 2018.01.15, 8 is the latest version.)
% 2 is the right one.
propSS = create_propCDR4basename('Acro','SS','Binning',propSP.binning,...
    'Wavelength_filter',propSP.wavelength_filter,'Side',propSP.side,...
    'SENSOR_ID',propSP.sensor_id,'Version',2);
% basenameSS = get_basenameCDR4_fromProp(propSS);
[dirpathSS,~,~,basenameSS] = get_dirpath_cdr_fromProp(propSS,varargin{:});
% basenameSSPtrn = get_basenameCDR4_fromProp(propSS);
% basenameSS = readDownloadBasename_v3(basenameSSPtrn,dirpathSS,remote_subdirSS);
if ~isempty(dirpathSS)
    % This will throw an error.
    basenameSS
    error('Please fix this');
    SSdata = CRISMdata(basenameSS,dirpathSS);
else
    error('%s cannot be found',basenameSS);
end

%VNIR SH data (as of 2018.01.15, 8 is the latest version.)
propSH = create_propCDR4basename('Acro','SH','Binning',propSP.binning,...
    'Wavelength_filter',propSP.wavelength_filter,'Side',propSP.side,...
    'SENSOR_ID',propSP.sensor_id,'Version',4);
[dirpathSH,~,~,basenameSH] = get_dirpath_cdr_fromProp(propSH);
% basenameSHPtrn = get_basenameCDR4_fromProp(propSH);
% basenameSH = readDownloadBasename_v3(basenameSHPtrn,dirpathSH,remote_subdirSH,varargin{:});
if ~isempty(dirpathSH)
    SHdata = CRISMdata(basenameSH,dirpathSH);
else
    error('%s cannot be found',basenameSH);
end

end
