function [index,prd] = cartesianProd(varargin)
%CARTPRODINDEX Cartesian product of sets
%   index = cartesianProd(sz) calculates the indices of the
%   cartesian product of arrays A_i with sz(i) = numel(A_i).
%   The final product can be calculated by:
%   [A_i(index(i,:))];
%
%   [index,otpt] = cartesianProd(A,B,...) automatically
%   finds the number of elements in each array and returns
%   the final product as a cell array
%
%   Cartesian product of any number of arrays is an array
%   which contains all possible tupples of elements in each
%   array. For example, The cartesian product of two sets A
%   = {1,2} and B = {a,b,c} is a 6 member set C =
%   {(1,a),(1,b),(1,3),(2,a),(2,b),(2,c)}.
%   In tensorial notation if A = B*C*D* ... then $A_{ijk...}
%   = B_i C_j D_k$
%
%   cartesianProd generates indices for final product and
%   then, if the input arrays are provided, it concats them
%   to produce the final product based on the generated
%   indices.
%
%   Example: Calculate the indices of cartesian product of
%   tree inputs with 2,3 and 3 number of elements
%   respectively.
%
%   index = cartesianProd([2,3,3])
%
%     index =
%
%       1     1     1
%       2     1     1
%       1     2     1
%       2     2     1
%       1     3     1
%       2     3     1
%       1     1     2
%       2     1     2
%       1     2     2
%       2     2     2
%       1     3     2
%       2     3     2
%       1     1     3
%       2     1     3
%       1     2     3
%       2     2     3
%       1     3     3
%       2     3     3
%
%   Example: Calculate the indices of cartesian product of
%   A=[1 2; 3 4] and B = ["a","b"]
%
%   A=[1 2; 3 4];
%   B = ["a","b"];
%   [index,prd] = cartesianProd([1 2; 3 4],["a","b"])
%
%       index =
%
%           1     1
%           2     1
%           3     1
%           4     1
%           1     2
%           2     2
%           3     2
%           4     2
%
%
%       prd =
%
%         8×2 cell array
%
%           {[1]}    {["a"]}
%           {[3]}    {["a"]}
%           {[2]}    {["a"]}
%           {[4]}    {["a"]}
%           {[1]}    {["b"]}
%           {[3]}    {["b"]}
%           {[2]}    {["b"]}
%           {[4]}    {["b"]}
%
%   Also the first and second colimns of prd are cell
%   equivalents of A and B rearranged by first and second
%   columns of index respectively:
%
%   A(index(:,1)) = [prd{:,1}]';
%
%
narginchk(1,inf);
if nargin == 1
    if ~isnumeric(varargin{1})
        error('cartesianProd:Input:Size:nonNumeric',...
            "The size of sets must be an array " + ...
            "of positive integer numbers");
    end
    sz = (varargin{1}(:)).';
    if any(sz<0 | sz ~= floor(sz))
        error('cartesianProd:Input:Size:nonPosInt',...
            "The size of sets must be an array " + ...
            "of positive integer numbers");
    end
else
    sz = cellfun(@numel,varargin);
    prd = cellfun(@(x)num2cell((x(:))'),varargin,...
        'uniformoutput',false);
end

index = repmat((0:prod(sz)-1)',1,numel(sz));
index = floor(index./[1,cumprod(sz(1:end-1))]);
index = mod(index,sz)+1;

if nargin == 1
    return;
end

prd = [prd{:}];
prd = prd(index + [0,cumsum(sz(1:end-1))]);


%% From License File
% 
% Copyright (c) 2021, saeed oveisi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% 
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
