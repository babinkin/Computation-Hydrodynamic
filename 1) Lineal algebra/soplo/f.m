function f = f(q, k, M0)
 f = q - M0*(((k+1)/2)^((k+1)/(2*(k -1)))) / ((1 + (k - 1) /2 * M0 ^ 2)^((k+1)/(2*(k -1))));
end