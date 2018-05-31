#!/usr/bin/ruby

stuff = <<STUFF;
1 1 1  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2 
1 1 2  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(%e^(r/2)*sinh(t/2))/(2*sqrt(r-1))+(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2+(%e^(r/2)*cosh(t/2))/(2*sqrt(r-1)) 
1 2 1  (sqrt(r-1)*%e^(r/2)*cosh(t/2))/2-(sqrt(r-1)*%e^(r/2)*sinh(t/2))/2 
1 2 2  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(%e^(r/2)*sinh(t/2))/(2*sqrt(r-1))-(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(r-1)) 
2 1 1  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2 
2 1 2  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2-(%e^(r/2)*sinh(t/2))/(2*sqrt(1-r))+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(1-r)) 
2 2 1  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2-(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2 
2 2 2  (-(sqrt(1-r)*%e^(r/2)*sinh(t/2))/2)+(%e^(r/2)*sinh(t/2))/(2*sqrt(1-r))+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(1-r)) 
3 1 1  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2 
3 1 2  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(%e^(r/2)*sinh(t/2))/(2*sqrt(r-1))+(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2+(%e^(r/2)*cosh(t/2))/(2*sqrt(r-1)) 
3 2 1  (sqrt(r-1)*%e^(r/2)*cosh(t/2))/2-(sqrt(r-1)*%e^(r/2)*sinh(t/2))/2 
3 2 2  (sqrt(r-1)*%e^(r/2)*sinh(t/2))/2+(%e^(r/2)*sinh(t/2))/(2*sqrt(r-1))-(sqrt(r-1)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(r-1)) 
4 1 1  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2 
4 1 2  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2-(%e^(r/2)*sinh(t/2))/(2*sqrt(1-r))+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(1-r)) 
4 2 1  (sqrt(1-r)*%e^(r/2)*sinh(t/2))/2-(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2 
4 2 2  (-(sqrt(1-r)*%e^(r/2)*sinh(t/2))/2)+(%e^(r/2)*sinh(t/2))/(2*sqrt(1-r))+(sqrt(1-r)*%e^(r/2)*cosh(t/2))/2-(%e^(r/2)*cosh(t/2))/(2*sqrt(1-r)) 
STUFF

stuff.each_line do |line|
  if line=~/(\d) (\d) (\d) (.*)/ then
    region,i,j,c = $1.to_i,$2.to_i,$3.to_i,$4
    c.gsub!("%e^","exp")
    c.gsub!("^","**")
    print "  # region #{region}\n    jacobian[#{i-1}][#{j-1}] = #{c}\n"
  end
end
