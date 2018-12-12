function gradf=fgradient_regular(f,order,d) 
% Compute gradient of a function on a regular grid with different order
% INPUTS
%   f    : 1d, 2d or 3d field
%   order: differentiation order, 2 for now
% Optional argument
%   d    : cell spacing in each dimensions, by default d = 1 or d=[1,1] or d=[1,1,1].
%
% AUTHOR: E.Branlard

%% TEST function
if nargin==0
    f=1:100;
    order=4;
    grad=fgradient_regular(f,order)
    return;
end

n    = setdiff(size(f),1); % Grid Dimensions
nDim = length(n);          % Physical dimension 1D, 2D, 3D

if ~exist('d','var'); d= ones(1,nDim); end;


gradf=zeros([nDim n(:)']);


d2=2*d;

if order==2 
    % --------------------------------------------------------------------------------
    % --- ORDER 2 
    % --------------------------------------------------------------------------------
    % left boundary: forward difference   (a=0, b=3)&
    % right boundary: backward difference (a=3, b=0)&
    % elsewhere: centered difference (a=1,b=1)&
    if nDim==1 
        for i = 1:n(1);
            % x-derivatives
            if (i==1) 
                df_dx     =  ( - 3 * f(i  ) + 4 * f(i+1) - 1 * f(i+2))/d2(1);
            elseif (i==n(1)) 
                df_dx     =  (   1 * f(i-2) - 4 * f(i-1) + 3 * f(i  ))/d2(1);
            else
                df_dx     =  ( - 1 * f(i-1) + 1 * f(i+1))/d2(1);
            end
            gradf(1,i)=df_dx;
        end % for
    elseif nDim==2 
        for j = 1:n(2);
        for i = 1:n(1);
            % x-derivatives
            if (i==1)
                df_dx     =  ( - 3 * f(i  ,j) + 4 * f(i+1,j) - 1 * f(i+2,j))/d2(1);
            elseif (i==n(1)) 
                df_dx     =  (   1 * f(i-2,j) - 4 * f(i-1,j) + 3 * f(i  ,j))/d2(1);
            else
                df_dx     =  ( - 1 * f(i-1,j) + 1 * f(i+1,j))/d2(1);
            end
            % y-derivatives
            if (j==1)
                df_dy     =  ( - 3 * f(i,j  ) + 4 * f(i,j+1)  - 1 * f(i,j+2))/d2(2);
            elseif (j==n(2)) 
                df_dy     =  (   1 * f(i,j-2) - 4 * f(i,j-1) + 3 * f(i,j  ))/d2(2);
            else
                df_dy     =  ( - 1 * f(i,j-1) + 1 * f(i,j+1))/d2(2);
            end
            gradf(1,i,j)=df_dx;
            gradf(2,i,j)=df_dy;
        end
        end % for
    elseif nDim==3 
        for k = 1:n(3);
        for j = 1:n(2);
        for i = 1:n(1);
            % x-derivatives
            if (i==1) 
                df_dx     =  ( - 3 * f(i  ,j,k) + 4 * f(i+1,j,k) - 1 * f(i+2,j,k))/d2(1);
            elseif (i==n(1)) 
                df_dx     =  (   1 * f(i-2,j,k) - 4 * f(i-1,j,k) + 3 * f(i  ,j,k))/d2(1);
            else
                df_dx     =  ( - 1 * f(i-1,j,k) + 1 * f(i+1,j,k))/d2(1);
            end
            % y-derivatives
            if (j==1) 
                df_dy     =  ( - 3 * f(i,j  ,k) + 4 * f(i,j+1,k)  - 1 * f(i,j+2,k))/d2(2);
            elseif (j==n(2)) 
                df_dy     =  (   1 * f(i,j-2,k) - 4 * f(i,j-1,k) + 3 * f(i,j  ,k))/d2(2);
            else
                df_dy     =  ( - 1 * f(i,j-1,k) + 1 * f(i,j+1,k))/d2(2);
            end
            % z-derivatives
            if (k==1) 
                df_dz     =  ( - 3 * f(i,j,k  ) + 4 * f(i,j,k+1) - 1 * f(i,j,k+2))/d2(3);
            elseif (k==n(3)) 
                df_dz     =  (   1 * f(i,j,k-2) - 4 * f(i,j,k-1) + 3 * f(i,j,k  ))/d2(3);
            else
                df_dz     =  ( - 1 * f(i,j,k-1) + 1 * f(i,j,k+1))/d2(3);
            end
            gradf(1,i,j,k)=df_dx;;
            gradf(2,i,j,k)=df_dy;;
            gradf(3,i,j,k)=df_dz;;
        end
        end
        end % for
    else
        error('Dimension not supported %d for order %d',nDim,order);
    end
elseif (order==4) 
    % --------------------------------------------------------------------------------
    % --- ORDER 4 
    % --------------------------------------------------------------------------------
    % left boundary: forward difference   (a=0, b=5) and (a=1, b=3)&
    % right boundary: backward difference (a=3, b=1) and (a=0, b=5)&
    % elsewhere: centered difference (a=2,b=2)&
    if nDim==1 
        for i = 1:n(1)
        % x-derivatives 
        if (i==1) 
            df_dx = ( - 25/6 * f( i)+ 8      * f( i+1)- 6    * f( i+2) +8/3 * f( i+3) -1/2   * f( i+4))/d2(1);
            elseif (i==2)
            df_dx= ( - 1/2   * f( i-1) - 5/3 * f( i) +3      * f( i+1)-1    * f( i+2)+ 1/6   * f( i+3))/d2(1);
            elseif (i==n(1)-1)
            df_dx= ( - 1/6   * f( i-3) +1    * f( i-2) -3    * f( i-1)+ 5/3 * f( i)+ 1/2     * f( i+1))/d2(1);
            elseif (i==n(1))
            df_dx = (1/2     * f( i-4)-8/3   * f( i-3)+ 6    * f( i-2) - 8  * f( i-1) + 25/6 * f( i))/d2(1);
            else
            df_dx= ( 1/6     * f( i-2) - 4/3 * f( i-1) + 4/3 * f( i+1)- 1/6 * f( i+2))/d2(1);
        end
        gradf(i)=df_dx;;
        end
    elseif nDim==2
        for j = 1:n(2)
        for i = 1:n(1)
        % x-derivatives 
        if (i==1) 
            df_dx = ( - 25/6 * f( i,j)+ 8      * f( i+1,j)- 6    * f( i+2,j) +8/3 * f( i+3,j) -1/2   * f( i+4,j))/d2(1);
            elseif (i==2)
            df_dx= ( - 1/2   * f( i-1,j) - 5/3 * f( i,j) +3      * f( i+1,j)-1    * f( i+2,j)+ 1/6   * f( i+3,j))/d2(1);
            elseif (i==n(1)-1)
            df_dx= ( - 1/6   * f( i-3,j) +1    * f( i-2,j) -3    * f( i-1,j)+ 5/3 * f( i,j)+ 1/2     * f( i+1,j))/d2(1);
            elseif (i==n(1))
            df_dx = (1/2     * f( i-4,j)-8/3   * f( i-3,j)+ 6    * f( i-2,j) - 8  * f( i-1,j) + 25/6 * f( i,j))/d2(1);
            else
            df_dx= ( 1/6     * f( i-2,j) - 4/3 * f( i-1,j) + 4/3 * f( i+1,j)- 1/6 * f( i+2,j))/d2(1);
        end
        % y-derivatives
        if (j==1)
            df_dy = ( - 25/6 * f( i,j)+ 8      * f( i,j+1)- 6    * f( i,j+2) +8/3 * f( i,j+3) -1/2   * f( i,j+4))/d2(2);
            elseif (j==2)
            df_dy= ( - 1/2   * f( i,j-1) - 5/3 * f( i,j) +3      * f( i,j+1)-1    * f( i,j+2)+ 1/6   * f( i,j+3))/d2(2);
            elseif (j==n(2)-1)
            df_dy= ( - 1/6   * f( i,j-3) +1    * f( i,j-2) -3    * f( i,j-1)+ 5/3 * f( i,j)+ 1/2     * f( i,j+1))/d2(2);
            elseif (j==n(2))
            df_dy = (1/2     * f( i,j-4)-8/3   * f( i,j-3)+ 6    * f( i,j-2) - 8  * f( i,j-1) + 25/6 * f( i,j))/d2(2);
            else
            df_dy= ( 1/6     * f( i,j-2) - 4/3 * f( i,j-1) + 4/3 * f( i,j+1)- 1/6 * f( i,j+2))/d2(2);
            end
        gradf(1,i,j)=df_dx;;
        gradf(2,i,j)=df_dy;;
        end
        end
    elseif nDim==3
        for k = 1:n(3)
        for j = 1:n(2)
        for i = 1:n(1)
        % x-derivatives 
        if (i==1) 
            df_dx = ( - 25/6 * f( i,j,k)+ 8 *      f( i+1,j,k)- 6 *    f( i+2,j,k) +8/3 *  f( i+3,j,k) -1/2 *   f( i+4,j,k))/d2(1);
        elseif (i==2)
            df_dx= ( - 1/2 *   f( i-1,j,k) - 5/3 * f( i,j,k) +3 *      f( i+1,j,k)-1 *     f( i+2,j,k)+ 1/6 *   f( i+3,j,k))/d2(1);
        elseif (i==n(1)-1)
            df_dx= ( - 1/6 *   f( i-3,j,k) +1 *    f( i-2,j,k) -3 *    f( i-1,j,k)+ 5/3 *  f( i,j,k)+ 1/2 *     f( i+1,j,k))/d2(1);
        elseif (i==n(1))
            df_dx = (1/2 *     f( i-4,j,k)-8/3 *   f( i-3,j,k)+ 6 *    f( i-2,j,k) - 8 *   f( i-1,j,k) + 25/6 * f( i,j,k))/d2(1);
        else
            df_dx= ( 1/6 *     f( i-2,j,k) - 4/3 * f( i-1,j,k) + 4/3 * f( i+1,j,k)- 1/6 *  f( i+2,j,k))/d2(1);
        end
        % y-derivatives
        if (j==1)
            df_dy = ( - 25/6 * f( i,j,k)+ 8 *      f( i,j+1,k)- 6 *    f( i,j+2,k) +8/3 *  f( i,j+3,k) -1/2 *   f( i,j+4,k))/d2(2);
        elseif (j==2)
            df_dy= ( - 1/2 *   f( i,j-1,k) - 5/3 * f( i,j,k) +3 *      f( i,j+1,k)-1 *     f( i,j+2,k)+ 1/6 *   f( i,j+3,k))/d2(2);
        elseif (j==n(2)-1)
            df_dy= ( - 1/6 *   f( i,j-3,k) +1 *    f( i,j-2,k) -3 *    f( i,j-1,k)+ 5/3 *  f( i,j,k)+ 1/2 *     f( i,j+1,k))/d2(2);
        elseif (j==n(2))
            df_dy = (1/2 *     f( i,j-4,k)-8/3 *   f( i,j-3,k)+ 6 *    f( i,j-2,k) - 8 *   f( i,j-1,k) + 25/6 * f( i,j,k))/d2(2);
        else
            df_dy= ( 1/6 *     f( i,j-2,k) - 4/3 * f( i,j-1,k) + 4/3 * f( i,j+1,k)- 1/6 *  f( i,j+2,k))/d2(2);
        end

        % z-derivatives
        if (k==1)
            df_dz = ( - 25/6 * f( i,j,k)+ 8 *      f( i,j,k+1)- 6 *    f( i,j,k+2) +8/3 *  f( i,j,k+3) -1/2 *   f( i,j,k+4))/d2(3);
        elseif (k==2)
            df_dz= ( - 1/2 *   f( i,j,k-1) - 5/3 * f( i,j,k) +3 *      f( i,j,k+1)-1 *     f( i,j,k+2)+ 1/6 *   f( i,j,k+3))/d2(3);
        elseif (k==n(3)-1)
            df_dz= ( - 1/6 *   f( i,j,k-3) +1 *    f( i,j,k-2) -3 *    f( i,j,k-1)+ 5/3 *  f( i,j,k)+ 1/2 *     f( i,j,k+1))/d2(3);

        elseif (k==n(3))
            df_dz= (1/2 *      f( i,j,k-4) -8/3 *  f( i,j,k-3) + 6 *   f( i,j,k-2)- 8 *    f( i,j,k-1)+ 25/6 *  f( i,j,k))/d2(3);
        else
            df_dz = ( 1/6 *    f( i,j,k-2)- 4/3 *  f( i,j,k-1)+ 4/3 *  f( i,j,k+1) - 1/6 * f( i,j,k+2))/d2(3);
        end
        gradf(1,i,j,k)=df_dx;;
        gradf(2,i,j,k)=df_dy;;
        gradf(3,i,j,k)=df_dz;;
        end
        end
        end % for
    else 
        error('nDim %d not implemented for order %d',nDim,order);
    end
else
   error('Order not implemented: %d',order);
end
