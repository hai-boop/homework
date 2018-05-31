function H= JacobianH(Xn)
    x=Xn(1);a =100;
    H=zeros(1,2);
    xa=1+(x/a)^2;
    H(1,1)=1/xa*1/a;H(1,2)=0;
