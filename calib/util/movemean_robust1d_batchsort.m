function [xhat] = movemean_robust1d_batchsort(x,varargin)
% [xhat] = movemean_robust1d(x)
%  sub routine for performing moving mean while removing the 4 extreme (2 
%  largest and 2 smallest) points. Currently window size is fixed at 15. 
%  
%  Note:
%  No padding is perfomed around edges, meaning the window size gets 
%  smaller (down to 8) around both of the edges. An exception occurs for 
%  the extreme points removal:
%  When the window size gets smaller than 10, only 2 extreme (the largest 
%  and the smallest) points are removed.
%  INPUT
%   x: 1-dimensional vector [Lx1]
%  OUTPUT
%   xhat: smoothened x (same size as x)
%  OPTIONAL PARAMETERS
%   'WINDOW_SIZE': window size of the moving mean
%       (default) 15

s = 15; % size of the window
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'WINDOW_SIZE'
                s = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error('Undefined option %s',varargin{i});
        end
    end
end
sh_b = floor((s-1)./2); % before
sh_a = floor(s./2); % after

L = length(x);
xhat = nan(size(x));

[VV,II] = sort(x,'ascend');
[~,IImap] = ismember(1:L,II);

% define the size of the window (smaller around edges)
idx_b = max((1:L)-sh_b,1); idx_a = min((1:L)+sh_a,L);
s_actual = abs(idx_a-idx_b)+1;

for i=1:L
%    if i==26
%        fprintf('%d\n',i);
%    end
   % define the number of outliers on one side; normally 2, but reduced
   % around the edge region.
   if s_actual(i) < 3
       n_outlier_one_side = 0;
   elseif s_actual(i) < 10
       n_outlier_one_side = 1;
   else
       n_outlier_one_side = 2;
   end
   % sorting the array, this implementation only takes 30% of time compared
   % with the latter formulation.
   % [v,ii] = sort(x(idx_b:idx_a),'ascend');
   idx = false(1,L);
   idx(IImap(idx_b(i):idx_a(i))) = true;
   % v = VV(IImap(idx_b(i):idx_a(i))); bug
   v = VV(idx);
   % exclude nans
   v_notnan = ~isnan(v); % ii_notnan = ii(v_notnan); 
   v_notnan = v(v_notnan); 
   v_rmean = mean(v_notnan((n_outlier_one_side+1):(end-n_outlier_one_side)));
   
   % another implementation for outlier removal this is much slower...
%    v_s = x(idx_b:idx_a);
%    [vmax,imax] = maxk(v_s,n_outlier_one_side);
%    [vmin,imin]  = mink(v_s,n_outlier_one_side);
%    idx_val = true(s_actual,1);
%    idx_val(imax) = false; idx_val(imin) = false;
%    v_rmean = nanmean(v_s(idx_val));
       
   
   
   xhat(i) = v_rmean;
   
end

end