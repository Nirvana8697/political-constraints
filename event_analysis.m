% Pre-default behavor of the model

% Expressing GDP policy under repayment as log deviation from simulation average

X = mean(mean(GDP(:,:)));
r = 1.01;


GDP_event = log(GDP)-log(mean(GDP_simul));

quarter = 40;

% Initial debt 

debt = 0.12;  %(.22 of average simulation GDP)

% bgrid index is

init = 167;
b_index = zeros(quarter+1,1);
M = zeros(quarter,1);
shock = zeros(quarter,1);
b_index(1) = init;
bgrid = double(bgrid);
GDP_arg = zeros(quarter,1);
tax_arg = zeros(quarter,1);
cons_arg = zeros(quarter,1);
Tr_arg = zeros(quarter,1);
gov_arg = zeros(quarter,1);
B_arg = zeros(quarter+1,1);
labor_arg = zeros(quarter,1);

% load data

data_spreads = xlsread('argentina_spreads',1);
data_rgdp = xlsread('argentina_gdp.xlsx',1);


Spreads_arg = zeros(1, quarter);
% Minimization


    for i = 1:quarter
    [M(i), shock(i)] = min(abs(GDP_event(:,b_index(i))- data_rgdp(i)));
    t = find(abs(bgrid-Bpol(shock(i),b_index(i)))<0.00000001);
    b_index(i+1) = t;
    GDP_arg(i) = GDP_event(shock(i),b_index(i));
    if (d(shock(i),b_index(i)) ~= 1) 
    Spreads_arg(i) = (1./q(shock(i),b_index(i+1)))^4-r;
    else 
        Spreads_arg(i) = 0.29;
    end 
    end
    
    for i = 1:quarter
        tax_arg(i) = tax(shock(i),b_index(i));
    end
   
    
    for i = 1:quarter
        B_arg(1) = -bgrid(b_index(1));
        if (d(shock(i),b_index(i))~=1)
            B_arg(i+1) = -Bpol(shock(i),b_index(i));
        else
            B_arg(i) = 0;
        end
            
    end
    
    for i = 1:quarter
        cons_arg(i) = cons(shock(i),b_index(i));
    end
    
    for i = 1:quarter
        labor_arg(i) = labor(shock(i),b_index(i));
    end
    
    for i = 1:quarter
        Tr_arg(i) = Tr(shock(i),b_index(i));
    end
    
    for i = 1:quarter
        gov_arg(i) = gov(shock(i),b_index(i));
    end
    
    

% GDP of Argentina 16 periods prior to crisis till crisis

% GDP_arg = zeros(quarter,1);
% 
% for i = 1:quarter
%     GDP_arg(i) = Z(shock(i))*labor_policy_new(shock(i),b_index(i));
% end


% figure(1)
plot(Spreads_arg)
  ylim([0,0.16])
  hold on
  plot(data_spreads)
%  figure(2)
%   plot(GDP_arg)
%   hold on
%   plot(data_rgdp)
%figure(3)
%plot(tax_arg)
% figure(4)
 %plot(cons_arg)
% figure(5)
 %plot(Tr_arg)
% figure(6)
%plot(gov_arg)
%plot(B_arg(1:quarter))
%plot(labor_arg)

