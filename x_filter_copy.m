function y = x_filter(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

y = zeros(size(x));
x1 = [0;x(1:end-1)];
y1 = 0;
for i = 1:length(x)
    y(i,1) = 0.9871*x(i,1) - 0.9871*x1(i,1) + 0.9742*y1;
    y1 = y(i,1);
end