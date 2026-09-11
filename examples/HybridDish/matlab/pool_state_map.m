function vp = pool_state_map(v, blockX, blockY)
%POOL_STATE_MAP  Average contiguous blocks on a 20x10 HybridDish map.
%   vp = POOL_STATE_MAP(v, blockX, blockY)
%   v is n-by-200. Voxel index matches Java: id = x*10 + y, x=0..19, y=0..9.
%   M(x+1, y+1) = v(x*10 + y + 1). No blur, no interpolation.
%
%   block 1x1 -> 200 cols (identity)
%   block 2x2 -> 50 cols  (10x5)
%   block 4x5 -> 10 cols  (5x2)
%   block 20x10 -> 1 col  (global mean)

    if nargin < 3
        error('pool_state_map:args', 'Need v, blockX, blockY');
    end
    nx = 20;
    ny = 10;
    if size(v, 2) ~= nx * ny
        error('pool_state_map:size', 'Expected %d map columns, got %d', nx*ny, size(v, 2));
    end
    if mod(nx, blockX) ~= 0 || mod(ny, blockY) ~= 0
        error('pool_state_map:block', 'Block %dx%d does not tile 20x10', blockX, blockY);
    end
    if blockX == 1 && blockY == 1
        vp = v;
        return
    end

    nOutX = nx / blockX;
    nOutY = ny / blockY;
    nRows = size(v, 1);
    vp = zeros(nRows, nOutX * nOutY);
    for r = 1:nRows
        M = reshape(v(r, :), [ny, nx])';  % 20 x 10, M(x+1,y+1)
        P = zeros(nOutX, nOutY);
        for i = 1:nOutX
            for j = 1:nOutY
                blk = M((i-1)*blockX+1:i*blockX, (j-1)*blockY+1:j*blockY);
                P(i, j) = mean(blk, 'all');
            end
        end
        vp(r, :) = reshape(P', 1, []);  % y-fast, x-outer
    end
end
