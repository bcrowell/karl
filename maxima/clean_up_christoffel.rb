#!/usr/bin/ruby

# This is to clean up maxima's output of the christoffel coefficients for K-S coordinates.
# There is output for python and latex. Can comment out one or the other as desired.
# The symbol om corresponds to what I'm now calling \ell, or r-1.
# Order of indices is ctensor's:
#    symmetric on 1st 2 indices
#    contravariant on final index

stuff = <<STUFF;
0 0 0 (sinh(a)*%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^((-b)+a-1))^2+(2\
*sinh(a)*%e^(lambert_w(%e^((-b)+a-1))+b+1)-%e^a*cosh(a))*lambert_w(%e^((-b)+a-\
1))+sinh(a)*%e^(lambert_w(%e^((-b)+a-1))+b+1)-2*%e^a*cosh(a))/(cosh(a)*%e^(lam\
bert_w(%e^((-b)+a-1))+b+1)+2*cosh(a)*%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert\
_w(%e^((-b)+a-1))+cosh(a)*%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^((-b)\
+a-1))^2) 
0 2 2 %e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+b\
+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^((\
-b)+a-1))^2) 
0 3 3 %e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+b\
+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^((\
-b)+a-1))^2) 
0 4 4 %e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+b\
+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^((\
-b)+a-1))^2) 
1 1 1 (%e^(lambert_w(%e^((-b)+a-1))+b+1)*sinh(b)*lambert_w(%e^((-b)+a-1))^2+(%\
e^a*cosh(b)+2*%e^(lambert_w(%e^((-b)+a-1))+b+1)*sinh(b))*lambert_w(%e^((-b)+a-\
1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*sinh(b)+2*%e^a*cosh(b))/(%e^(lambert_w(%\
e^((-b)+a-1))+b+1)*cosh(b)+2*%e^(lambert_w(%e^((-b)+a-1))+b+1)*cosh(b)*lambert\
_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*cosh(b)*lambert_w(%e^((-b)\
+a-1))^2) 
1 2 2 -%e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+\
b+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^(\
(-b)+a-1))^2) 
1 3 3 -%e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+\
b+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^(\
(-b)+a-1))^2) 
1 4 4 -%e^a/(%e^(lambert_w(%e^((-b)+a-1))+b+1)+2*%e^(lambert_w(%e^((-b)+a-1))+\
b+1)*lambert_w(%e^((-b)+a-1))+%e^(lambert_w(%e^((-b)+a-1))+b+1)*lambert_w(%e^(\
(-b)+a-1))^2) 
2 2 0 -(%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
2 2 1 (%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
3 3 0 -(%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
3 3 1 (%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
4 4 0 -(%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
4 4 1 (%e^-b*(%e^a*lambert_w(%e^((-b)+a-1))+%e^a))/(2*cosh(a)*cosh(b)) 
STUFF

###################################################################

stuff.gsub!(/\\\n/,'') # merge broken lines into single long lines

###################################################################
# python output
###################################################################
coord = ['a','b','i','j','k']

stuff.each_line do |line|
  if line=~/(\d) (\d) (\d) (.*)/ then
    i,j,k,c = $1.to_i,$2.to_i,$3.to_i,$4
    c.gsub!("%e^a","exp(a)")
    c.gsub!("%e^b","exp(b)")
    c.gsub!("%e^-a","exp(-a)")
    c.gsub!("%e^-b","exp(-b)")
    c.gsub!("%e^","exp")
    c.gsub!("^","**")
    if i!=j then top=2 else top=1 end
    1.upto(top) { |m|
      if m==2 then u=i; i=j; j=u end
      print "  ch[#{i}][#{j}][#{k}] = #{c}\n  #   ... ^#{coord[k]} _#{coord[i]} #{coord[j]}\n"
    }
  end
end


###################################################################
# latex output
###################################################################

if true then

coord = ['a','b','i','j','k']

stuff.each_line do |line|
  if line=~/(\d) (\d) (\d) (.*)/ then
    i,j,k,c = $1.to_i,$2.to_i,$3.to_i,$4
    c.gsub!("%e^a","e^a ")
    c.gsub!("%e^b","a^b ")
    c.gsub!("%e^-a","e^{-a}")
    c.gsub!("%e^-b","e^{-b}")
    c.gsub!("%e^","exp")
    #c.gsub!("lambert_w(exp((-b)+a-1))","W(e^{-b+a-1})")
    c.gsub!("lambert_w(exp((-b)+a-1))","Q")
    c.gsub!("exp(Q+b+1)","P")
    c.gsub!(/(exp|cosh|sinh)/) {"\\#{$1}"}
    c.gsub!("\*","")
    c.gsub!("%e^","e^")
    c.gsub!("lambert_w","W")
    print "\\Gamma\\indices{^#{coord[k]}_#{coord[i]}_#{coord[j]}} &= #{c} \\\\\n"
  end
end

end
