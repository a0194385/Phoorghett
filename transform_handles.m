@(x) conv(x(:),[1 -1],'same'); % first line is transform and second line is inverse transform.
@(x) flipud(-cumsum(flipud(x(:))));