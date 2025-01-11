function y = quantization(x, bit, max)
% 'max' is supposed be max(abs(signal)).
% Input signal 'x' must be less than or equal to 'max'. Otherwise, quantization outputs with more bits than 'bit' will be generated.

if bit == inf
    y = x;
else
    qStep = max/2^(bit-1);
    y_real = floor(real(x)/qStep)*qStep + qStep/2;
    y_imag = floor(imag(x)/qStep)*qStep + qStep/2;
    y = y_real + 1j*y_imag;
end
disp('')
