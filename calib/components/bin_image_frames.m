function [img_fr_bn] = bin_image_frames(img_fr,varargin)
% [img_fr_bn] = bin_image_frames(img_fr,varargin)
%  apply binning to image_frames 
%   Input
%     img_fr   : image frames to be binned [L x S x B]
%   Output Parameters
%     img_fr_bn: binned image frames [L x Sb x B] (Sb = S/binx)
%   Optional Parameters
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
%    If you only specify 'BINNING', 'BINX' is obtained automatically. If
%    you specify 'BINX', whatever 'BINNING' parameter will not be
%    considered.
% [~,S,Band] = size(RSPj);

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end


S_b = crism_getSampleSize(binx);
[L,S_ori,B] = size(img_fr);

if S_b==S_ori
    img_fr_bn = img_fr;
else
    h = fspecial('average',[binx,1]); % convolution vector
    img_fr_bn = nan([L,S_b,B]);
    for l=1:L
        image_frames_l = squeeze(img_fr(l,:,:)); %[S x B]
        image_frames_lbn = conv2(image_frames_l,h,'valid'); % take convolution and 
        image_frames_lbn = image_frames_lbn(1:binx:end,:); % downsample
        img_fr_bn(l,:,:) = reshape(image_frames_lbn,[1,S_b,B]);
    end

end

end