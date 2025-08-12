function pos = axgridpos(nf,mf,i,pdx,pdy,varargin)
% Function to calculate position of axes in a grid given dimensions of grid
% and inner/outer padding values. This is meant to give better results than
% the native "subplot." The function returns a position vector which can
% then be used to create axes.
%
% pos = axgridpos(nf,mf,i,pdx,pdy)
% pos = axgridpos(nf,mf,i,pdx,pdy,pdx_outer,pdy_outer)
% pos = axgridpos(nf,mf,i,pdx,pdy,_,_,flip)
%
% Input
%   nf: number of rows
%   mf: number of columns
%   i: index of axes (defaults to row-wise counting)
%   pdx: x padding in normalized units
%   pdy: y padding in normalized units
%   pdx_outer: (optional) x padding between grid and edge of figure
%   pdy_outer: (optional) y padding between grid and edge of figure
%   flip: (optional) input 'flip' to count axes column-wise
%
% KJW
% 17 Jan 2023

% catch some errors
if (round(nf)~=nf || nf<1) || (round(mf)~=mf || mf<1)
    error('"nf" and "mf" must be positive integers')
end

if round(i)~=i || i<1 || i>nf*mf
    error('"i" must be an integer between 1 and nf*mf')
end

% parse varargin
pdx_outer = pdx;
pdy_outer = pdy;
if nargin==7 || nargin==8
    pdx_outer = varargin{1};
    pdy_outer = varargin{2};
end

flip = false;
if any(strcmp(varargin,'flip'))
    flip = true;
end

% coordinates of axes in the grid (bottom left is (0,0))
if ~flip
    p = mod((i-1),mf); % horizontal
    q = floor((mf*nf-i)/mf); % vertical
else % flip
    p = floor((i-1)/nf);
    q = mod(mf*nf-i,nf);
end

% height and width
width = (1-2*pdx_outer-(mf-1)*pdx)/mf;
height = (1-2*pdy_outer-(nf-1)*pdy)/nf;

% lower left position
x = pdx_outer + p*(width+pdx);
y = pdy_outer + q*(height+pdy);

pos = [x y width height];

