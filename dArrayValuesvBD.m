% Order for dArray:order 1:D11(Jacobian)D20 (related to variance) order 2: D12-D21-D30 order 3:D22-D13-D40-D31
%M

function dArray=dArrayValuesvBD(phi)

dArray=cell(1,9);

dArray{1}=phi/(phi-1); dArray{2}=2*phi; dArray{3}=-2/(1-phi); 
dArray{4}=(2-3*phi)/(1-phi); dArray{5}=0; dArray{6}=-2/(1-phi);
dArray{7}=0; dArray{8}=2*phi; dArray{9}=phi/(phi-1);   

end




