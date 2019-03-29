syms h p H r a

assume(h > 0)
assume(h < 1)
assume(p > 1);
assume(r > 1);
assume(a > 0);
assume(h - h^(p+1/2) > 0)

I = 1 / (2*h^(p+1/2)) * int(H^r, H, h-h^(p+1/2), h+h^(p+1/2));

I2 = 1 / (2*h^(p+1/2)) * int((H^r - h^r)^2, H, h-h^(p+1/2), h+h^(p+1/2));

I3 = 1 / (2*h^(p+1/2)) * int(H*exp(-a/H), H, h-h^(p+1/2), h+h^(p+1/2));