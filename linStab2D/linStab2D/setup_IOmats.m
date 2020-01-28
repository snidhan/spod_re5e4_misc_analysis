function [B,C]   = setup_IOmats(grid,io)

B_win       = reshape(sparse(grid.x>io.input_x_min&grid.x<io.input_x_max&grid.r>io.input_r_min&grid.r<io.input_r_max),grid.NrNx,1);

B = [ ...
    diag(B_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                diag(B_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(B_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(B_win).*speye(grid.NrNx)     0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(B_win).*speye(grid.NrNx)    
    ];

C_win       = reshape(sparse(grid.x>io.output_x_min&grid.x<io.output_x_max&grid.r>io.output_r_min&grid.r<io.output_r_max),grid.NrNx,1);

C = [ ...
    diag(C_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                diag(C_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(C_win).*speye(grid.NrNx)     0*speye(grid.NrNx)                0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(C_win).*speye(grid.NrNx)     0*speye(grid.NrNx)
    0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                0*speye(grid.NrNx)                diag(C_win).*speye(grid.NrNx)    
    ];