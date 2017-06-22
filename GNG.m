function [grn_val, test_vec] = GNG
% Gaussina Random Number Generator
% variables in capital are bit vectors
% variables in small letters are double scalars

% Uniform Random Number Generator
URNG_OUT = double(rand(1,64)>0.5);

% Bit assignment
SIGN = URNG_OUT(1);
OFFSET = URNG_OUT(2:3);
LZD_INPUT = URNG_OUT(4:64);

% Leading Zero Detector (5 bits)
% 61 with the probability of 1/2
% 60 with the probability of 1/4
% .
% .
% 30 with the probability of 1/2^32
lz_out = find(LZD_INPUT==0);
if isempty(lz_out) || lz_out(end)<30
    lz_out = 30;
else
    lz_out = lz_out(end);
end
%
seg_idx = 4*(61-lz_out+1)+vec2dec(OFFSET); % Integer 1~128
SEGMENT = de2bi(seg_idx,7); % 7 bit vectors

MASK_IN = URNG_OUT(4:18);
MASK_OUT = MASK_IN; MASK_OUT(lz_out:end)=0;
MULT_IN = fliplr(MASK_OUT); % QF = 15

% 128 segment range table generation.
% 1st column is the start and 2nd point is the end point of the segment.
A = cumsum(0.5./(2.^(0:32))); A = repmat(A,4,1); 
B = 0.25./2.^(0:32)/4;
table = A(:)+kron(B,[0:3]).';
segment_table = [table(1:128) table(2:129)];

% coefficient generation
% polynomical approximation to ICDF within each segment
coeff_table = zeros(128,3);
x = (0:1/2^15:1).';
for ii = 1:size(coeff_table,1)
    coeff_table(ii,:) =  pinv([x.^2 x ones(size(x))])*...
        norminv(segment_table(ii,1)+x*diff(segment_table(ii,:)));
end

% Polynomial Interpolation
seg_len = diff(segment_table(seg_idx,:));
% grn_val = polyval(coeff_table(seg_idx,:),seg_len/2^15*vec2dec(MULT_IN));
fixpt_coef = round(coeff_table(seg_idx,:)*2^18);
COEF0 = de2bi(fixpt_coef(3),21); % QF=18
COEF1 = de2bi(fixpt_coef(2),18); % QF=18
COEF2 = de2bi(fixpt_coef(1),18); % QF=18
interm_res = seg_len/2^15*vec2dec(MULT_IN)/2^15*vec2dec(COEF2)/2^18+vec2dec(COEF1)/2^18;
INTERM = de2bi(round(interm_res*2^18),18); % QF=18
poly_out = seg_len/2^15*vec2dec(MULT_IN)/2^15*vec2dec(INTERM)/2^18+vec2dec(COEF0)/2^18;
POLYOUT = de2bi(round(poly_out*2^13),16); % QF=13
if SIGN==1
    grn_val = -poly_out;
    GRNOUT = de2bi(vec2dec(1-POLYOUT)+1); % 2's complement
else
    grn_val = poly_out;
    GRNOUT = POLYOUT;
end
test_vec.URNG_OUT = URNG_OUT;
test_vec.SIGN = SIGN;
test_vec.OFFSET = OFFSET;
test_vec.LZD_INPUT = LZD_INPUT;
test_vec.SEGMENT = SEGMENT;
test_vec.MASK_IN = MASK_IN;
test_vec.MASK_OUT = MASK_OUT;
test_vec.MULT_IN = MULT_IN;
test_vec.SEGMENT = SEGMENT;
test_vec.COEF0 = COEF0;
test_vec.COEF1 = COEF1;
test_vec.COEF2 = COEF2;
test_vec.INTERM = INTERM;
test_vec.POLYOUT = POLYOUT;
test_vec.GRNOUT = GRNOUT;

function y = vec2dec(v)
y = sum(2.^(0:length(v)-1).*v);