function nout = normal_rand(mu,sig,minv,maxv,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         in - input 1
%         in2 - input 2
%
% OUTPUT:
%    p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
% deal with input variables

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    nout = [];
    ncount = 0;
    while ncount<num
        val = normrnd(mu,sig,1,1);
        if val>minv && val<maxv
            nout = [nout; val];
            ncount = ncount + 1;
        end
    end
    nout = nout(:);

    
    



































