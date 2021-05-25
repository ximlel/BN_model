syms h h_R q g Z_L Z_R;
fun = 2*(h-h_R)*(1-q^2*h/g/(2*h-h_R)^2/h_R^2)+Z_L-Z_R;
diff(fun,h_R)