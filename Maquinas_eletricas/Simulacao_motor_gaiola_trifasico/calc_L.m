function [Ls , Lm , Lr , L0] = calc_L(Xls , Xlr , Xm, Wsync)
Ls = (Xls + Xm) / Wsync;
Lm = Xm / Wsync;
Lr = (Xlr + Xm) / Wsync;

L0 = sqrt((Ls*Lr)- (Lm^2));

end

