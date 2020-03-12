function varargout = gradient(f,varargin)
%GRADIENT Approximate gradient.
%   [FX,FY] = GRADIENT(F) returns the numerical gradient of the
%   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) 
%   direction. FY corresponds to dF/dy, the differences in y (vertical) 
%   direction. The spacing between points in each direction is assumed to 
%   be one. When F is a vector, DF = GRADIENT(F) is the 1-D gradient.
%
%   [FX,FY] = GRADIENT(F,H), where H is a scalar, uses H as the
%   spacing between points in each direction.
%
%   [FX,FY] = GRADIENT(F,HX,HY), when F is 2-D, uses the spacing
%   specified by HX and HY. HX and HY can either be scalars to specify
%   the spacing between coordinates or vectors to specify the
%   coordinates of the points.  If HX and HY are vectors, their length
%   must match the corresponding dimension of F.
%
%   [FX,FY,FZ] = GRADIENT(F), when F is a 3-D array, returns the
%   numerical gradient of F. FZ corresponds to dF/dz, the differences
%   in the z direction. GRADIENT(F,H), where H is a scalar, 
%   uses H as the spacing between points in each direction.
%
%   [FX,FY,FZ] = GRADIENT(F,HX,HY,HZ) uses the spacing given by
%   HX, HY, HZ. 
%
%   [FX,FY,FZ,...] = GRADIENT(F,...) extends similarly when F is N-D
%   and must be invoked with N outputs and either 2 or N+1 inputs.
%
%   Note: The first output FX is always the gradient along the 2nd
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the 1st dimension of F, going across rows.  For the
%   third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%
%   Examples:
%       [x,y] = meshgrid(-2:.2:2, -2:.2:2);
%       z = x .* exp(-x.^2 - y.^2);
%       [px,py] = gradient(z,.2,.2);
%       contour(z), hold on, quiver(px,py), hold off
%
%   Class support for input F:
%      float: double, single
%
%   See also DIFF, DEL2.

%   Copyright 1984-2013 The MathWorks, Inc.

[f,ndim,loc,rflag] = parse_inputs(f,varargin);
nargoutchk(0,ndim);

% Loop over each dimension. 

varargout = cell(1,ndim);
siz = size(f);
% first dimension 
g  = zeros(size(f),class(f)); % case of singleton dimension
h = loc{1}(:); 
n = siz(1);
% Take forward differences on left and right edges
if n > 2
   g(1,:) = (-1/2*f(3,:)+2*f(2,:) - 3/2*f(1,:))/(h(2)-h(1));
   g(2,:) = (-1/2*f(4,:)+2*f(3,:) - 3/2*f(2,:))/(h(2)-h(1));
   g(n,:) = (3/2*f(n,:) - 2*f(n-1,:) +1/2*f(n-2,:))/(h(end)-h(end-1));
   g(n-1,:) = (3/2*f(n-1,:) - 2*f(n-2,:) +1/2*f(n-3,:))/(h(end)-h(end-1));
end

% Take centered differences on interior points
if n > 4
   h = h(5:n) - h(1:n-4);
   g(3:n-2,:) = bsxfun(@rdivide,(-1/6*f(5:n,:)+4/3*f(4:n-1,:)-4/3*f(2:n-3,:)+1/6*f(1:n-4,:)),h);
end

varargout{1} = g;

% second dimensions and beyond
for k = 2:ndim
   n = siz(k);
   newsiz = [prod(siz(1:k-1)) siz(k) prod(siz(k+1:end))];
   nf = reshape(f,newsiz);
   h = loc{k}(:).';   
   g  = zeros(size(nf),class(nf)); % case of singleton dimension

   % Take forward differences on left and right edges
   if n > 2
      g(:,1,:) = (-1/2*nf(:,3,:) + 2*nf(:,2,:) - 3/2*nf(:,1,:))/(h(2)-h(1));
      g(:,n,:) = (3/2*nf(:,n,:) -2*nf(:,n-1,:) + 1/2*nf(:,n-2,:))/(h(end)-h(end-1));
      g(:,2,:) = (-1/2*nf(:,4,:) + 2*nf(:,3,:) - 3/2*nf(:,2,:))/(h(2)-h(1));
      g(:,n-1,:) = (3/2*nf(:,n-1,:) -2*nf(:,n-2,:) + 1/2*nf(:,n-3,:))/(h(end)-h(end-1));
   end

   % Take centered differences on interior points
   if n > 4
      h = h(5:n) - h(1:n-4);
      g(:,3:n-2,:) = bsxfun(@rdivide,(-1/6*nf(:,5:n,:)+4/3*nf(:,4:n-1,:)-4/3*nf(:,2:n-3,:)+1/6*nf(:,1:n-4,:)),h);
   end

   varargout{k} = reshape(g,siz);
end 

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim > 1
    varargout(2:-1:1) = varargout(1:2);
elseif rflag
    varargout{1} = varargout{1}.';
end


%-------------------------------------------------------
function [f,ndim,loc,rflag] = parse_inputs(f,v)
%PARSE_INPUTS
%   [ERR,F,LOC,RFLAG] = PARSE_INPUTS(F,V) returns the spacing
%   LOC along the x,y,z,... directions and a row vector
%   flag RFLAG. 

loc = {};

% Flag vector case and row vector case.
ndim = ndims(f);
vflag = false;
rflag = false;
if isvector(f)
    ndim = 1;
    vflag = true;
    if isrow(f) % Treat row vector as a column vector
        rflag = true;
        f = f.';
    end
end;

indx = size(f);

% Default step sizes: hx = hy = hz = 1
if isempty(v)
    % gradient(f)
    loc = cell(1, ndims(f));
    for k = 1:ndims(f)
        loc(k) = {1:indx(k)};
    end
elseif isscalar(v) % gradient(f,h)
    % Expand scalar step size
    if isscalar(v{1})
        loc = cell(1, ndims(f));
        for k = 1:ndims(f)
            h = v{1};
            loc(k) = {h*(1:indx(k))};
        end
        % Check for vector case
    elseif vflag
        loc(1) = v(1);
    else
        error(message('MATLAB:gradient:InvalidInputs'));
    end
elseif ndims(f) == numel(v)  % gradient(f,hx,hy,hz,...)
    % Swap 1 and 2 since x is the second dimension and y is the first.
    loc = v;
    if ndim > 1
        loc(2:-1:1) = loc(1:2);
    end
    % replace any scalar step-size with corresponding position vector
    for k = 1:ndims(f)
        if isscalar(loc{k})
            loc{k} = loc{k}*(1:indx(k));
        end
    end 
else
    error(message('MATLAB:gradient:InvalidInputs'));
end
