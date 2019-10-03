
function pic2mov(n)

    clc; clf;
    vidObj = VideoWriter('ackley.avi');
    vidObj.Quality = 95;
    vidObj.FrameRate = 1;
    open(vidObj);

    % Create movie.
    for i = 0:1:n
        fn = sprintf('ackley/ackley_CHC_%d', i);
        data = imread([fn, '.png']);
        j = i+1;
        writeVideo(vidObj,data);
    end

    % Create AVI file.
    close(vidObj);

end

