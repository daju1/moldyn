function [km] = korrelation_momentum(i, v1, v2)

s1 = length(v1)
s2 = length(v2)

s = min(s1,s2)

km = 0;
k = 0;

for j = 1 : s - i
    km = km + v1(j+i) * v2(j);
    
    k = k + 1;    
end

km = km / k;
