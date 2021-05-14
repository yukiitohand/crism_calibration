function [Xhat] = crmcal_movemean_robust(X,dim)
% [Xhat] = crmcal_movemean_robust(X,dim)
%  Performing moving mean while removing the 4 extreme (2 largest and 2 
%  smallest) points. Currently window size is fixed at 15. 
%  
%  Note:
%  No padding is perfomed around edges, meaning the window size gets 
%  smaller (down to 8) around both of the edges. An exception occurs for 
%  the extreme points removal:
%  When the window size gets smaller than 10, only 2 extreme (the largest 
%  and the smallest) points are removed.
%
%  INPUT
%   X: 3-dimensional matrix [LxSxB]
%   dim: dimesion along which moving mean is performed. currently dim=1 is
%   only implemented.
%  OUTPUT
%   xhat: smoothened x (same size as x)
[~,~,B] = size(X);
Xhat = X;
switch dim
    case 1
        parfor bi=1:B
            Xhat(:,:,bi) =  movemean_robust2d(X(:,:,bi));
        end
    otherwise
        error('Not implemented');
end

end

function [Xhat_bi] = movemean_robust2d(X_bi)            
    [~,S] = size(X_bi);
    Xhat_bi = nan(size(X_bi));
    for si=1:S
        % [VV,II] = sort(X_bi,1,'ascend');
        Xhat_bi(:,si) = crmcal_movemean_robust1d_batchsort(X_bi(:,si));% ,VV(:,si),II(:,si));       
    end
end

function [xhat] = movemean_robust1d(x)
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
s = 15; % size of the window
sh_b = floor((s-1)/2); % before
sh_a = floor(s/2); % after

L = length(x);
xhat = [];
for i=1:L
   % define the size of the window (smaller around edges)
   idx_b = max(1,i-sh_b); idx_a = min(L,i+sh_a);
   s_actual = idx_b-idx_a+1;
   % define the number of outliers on one side; normally 2, but reduced
   % around the edge region.
   if s_actual < 10
       n_outlier_one_side = 1;
   else
       n_outlier_one_side = 2;
   end
   
   % sorting the array, this implementation only takes 30% of time compared
   % with the latter formulation.
   [v,ii] = sort(x(idx_b:idx_a),'ascend');
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