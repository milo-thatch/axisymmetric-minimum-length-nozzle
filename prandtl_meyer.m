%Ludovico Fossà 01/2021

function nu=prandtl_meyer(mach,gamma)
    nu=sqrt((gamma+1)./(gamma-1)).*atan(sqrt((gamma-1)./(gamma+1).*(mach.^2-1)))-atan(sqrt(mach.^2-1));
end