% This script is to load the policy functions for the linear model and test
% the accuracy of the results.

noStates = 201;
bsiz = 100;


Bpol10 = reshape(dlmread('Bpol_linear.txt'),[noStates, bsiz]);

tax10 = reshape(dlmread('tax_linear.txt'),[noStates, bsiz]);
 
q10 = reshape(dlmread('q_linear.txt'),[noStates, bsiz]);
 
Tr10 = reshape(dlmread('Tr_linear.txt'),[noStates, bsiz]);
 
d10 = reshape(dlmread('d_linear.txt'),[noStates, bsiz]);

% srep = reshape(dlmread('srep_linear.txt'),[noStates, bsiz, m]);
% 
% sdef = dlmread('sdef_linear.txt');
% 
% Proprep = reshape(dlmread('Proprep_linear.txt'),[noStates, bsiz]);
% 
% Propdef = reshape(dlmread('Propdef_linear.txt'),[noStates, bsiz]);
% 
% Prop = reshape(dlmread('Prop_linear.txt'),[noStates, bsiz]);
 
 gov10 = reshape(dlmread('g_linear.txt'),[noStates, bsiz]);
 
 GDP10 = reshape(dlmread('GDP_linear.txt'),[noStates, bsiz]);
 
 cons10 = reshape(dlmread('cons_linear.txt'),[noStates, bsiz]);
 
 labor10 = reshape(dlmread('labor_linear.txt'),[noStates, bsiz]);
 
 Vc10 = reshape(dlmread('Vc_linear.txt'),[noStates, bsiz]);
 
 Vd10 = dlmread('Vd_linear.txt');
 
 Bdef10 = dlmread('Bdef_linear.txt');
 
 govdef10 = dlmread('govdef_linear.txt');
 
 taxdef10 = dlmread('taxdef_linear.txt');
 
 
 
  %Vhyp = Vd1 - Bdef1;
 
 % Plots 
 
%   plot(-bgrid, Vc5(50,:))
%   ylim([-45,-25]);
%   hold on
%   yline(Vd5(50));
%   hold on
%   plot(-bgrid, Vc1(50,:))
%   hold on
%   yline(Vd1(50));
%   hold on
%   yline(Vhyp(50));
%   hold off
%  
%  
 

 
