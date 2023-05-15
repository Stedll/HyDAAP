% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   12/2015
%   contact: stefan.haerer@boku.ac.at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str_out] = Unix_Slash_PRACTISE(str_in)
%   Name:       Unix_Slash_PRACTISE
%   Purpose:    Check if PRACTISE is run in Unix environment and change 
%               in case the folder separator from "\" to "/"
%   
%
%   Input:      str_in (Input string with "\")
%   Output:     str_out (Output string with "/")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
persistent CacheVal;  % speeds up repeated calls
%
if isempty (CacheVal)
    CacheVal = (isunix () > 0);
end
%
if CacheVal==1
    str_out=strrep(str_in,'\','/');
else
    str_out=str_in;
end