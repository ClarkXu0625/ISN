function output = export_bistable(N, Nvec1, Nvec2, outputmat)
    output = zeros(Nvec1*(N+1) , Nvec2);
    for i = 0:N
        i_start = i*Nvec1 + 1;
        i_end = (i+1)*Nvec1;
        output(i_start:i_end, : ) = reshape(outputmat(i+1,:,:), [Nvec1, Nvec2]);
    end
end