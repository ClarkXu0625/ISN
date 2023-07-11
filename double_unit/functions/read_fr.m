function fr_vec = read_fr(N, Nt, output, i, j)
% read fr_vec(N x Nt) from an output vector, which can be further
% accessed with i (1<=i<=Nvec1) and j (1<=j<=Nvec2)
% Output vector is an 2-D vector that can be stored in a separate csv file
    i_start = N*(i-1)+1;
    i_end = N*i;
    j_start = Nt*(j-1)+1;
    j_end = Nt*j;

    fr_vec = output(i_start:i_end, j_start:j_end);