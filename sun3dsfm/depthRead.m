function depth = depthRead(filename, bit_shift)
    depth = imread(filename);
    if bit_shift
        depth = bitor(bitshift(depth,-3), bitshift(depth,16-3));
    end
    depth = single(depth)/1000;
end
