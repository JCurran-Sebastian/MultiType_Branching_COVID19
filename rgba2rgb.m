function [RGB_colour] = rgba2rgb(RGB_background,RGBA_color)
%rgba2rgb Get rgb colour map that looks similar to base colour with
%transparency applied

    % Extract alpha value from input
    alpha = RGBA_color(4);
    RGB_colour = zeros(1,3);
    for idx_itr = 1:length(RGB_colour)
        RGB_colour(idx_itr) = (1 - alpha) * RGB_background(idx_itr) + alpha * RGBA_color(idx_itr);
    end
end

% Based on following code demo: http://marcodiiga.github.io/rgba-to-rgb-conversion
% function rgba2rgb(RGB_background, RGBA_color)
% {
%     var alpha = RGBA_color.a;
% 
%     return new Color(
%         (1 - alpha) * RGB_background.r + alpha * RGBA_color.r,
%         (1 - alpha) * RGB_background.g + alpha * RGBA_color.g,
%         (1 - alpha) * RGB_background.b + alpha * RGBA_color.b
%     );
% }