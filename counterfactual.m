% Counterfactual GDP

GDP_event = log(GDP10)-log(mean(GDP_simul10));

B_arg10 = zeros(quarter+1,1);
tax_arg10 = zeros(quarter,1);
cons_arg10 = zeros(quarter,1);
labor_arg10 = zeros(quarter,1);
gov_arg10 = zeros(quarter,1);
Tr_arg10 = zeros(quarter,1);
Spreads_arg10 = zeros(quarter,1);
b_index = zeros(quarter+1,1);
init = 167;
b_index(1) = init;
GDP_arg10 = zeros(quarter,1);

for i = 1:quarter
        t = find(abs(bgrid-Bpol10(shock(i),b_index(i)))<0.00000001);
        b_index(i+1) = t;
end

for i = 1:quarter
        GDP_arg10(i) = GDP_event(shock(i),b_index(i));
end

for i = 1:quarter
        tax_arg10(i) = tax10(shock(i),b_index(i));
end

for i = 1:quarter
        cons_arg10(i) = cons10(shock(i),b_index(i));
end

for i = 1:quarter
        labor_arg10(i) = labor10(shock(i),b_index(i));
end

for i = 1:quarter
        gov_arg10(i) = gov10(shock(i),b_index(i));
end

for i = 1:quarter
        Tr_arg10(i) = Tr10(shock(i),b_index(i));
end

for i = 1:quarter
        B_arg10(1) = -bgrid(b_index(1));
        if (d10(shock(i),b_index(i))~=1)
            B_arg10(i+1) = -Bpol10(shock(i),b_index(i));
        else
            B_arg10(i) = 0;
        end
            
end

for i = 1:quarter
        if (d10(shock(i),b_index(i)) ~= 1) 
        Spreads_arg10(i) = (1./q10(shock(i),b_index(i+1)))^4-r;
        else 
        Spreads_arg10(i) = 0.29;
        end 
            
end


%plot(GDP_arg10)
plot(Spreads_arg10(1:quarter))