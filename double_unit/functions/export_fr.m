function output = export_fr(N, Nt, input_vec, output, i, j)
% store input_vec(N x Nt) into an output vector, which can be further
% accessed with i (1<=i<=Nvec1) and j (1<=j<=Nvec2)
% Output vector is an 2-D vector that can be stored in a separate csv file
    i_start = N*(i-1)+1;
    i_end = N*i;
    j_start = Nt*(j-1)+1;
    j_end = Nt*j;
    
    output(i_start:i_end, j_start:j_end) = input_vec;