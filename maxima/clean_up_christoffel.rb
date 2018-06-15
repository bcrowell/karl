#!/usr/bin/ruby

# This is to clean up maxima's output of the christoffel coefficients for K-S coordinates.
# There is output for python and latex. Can comment out one or the other as desired.
# The symbol om corresponds to what I'm now calling \ell, or r-1.
# Order of indices is ctensor's:
#    symmetric on 1st 2 indices
#    contravariant on final index

stuff = <<STUFF;
0 0 0 (sinh(a)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(\
a)*sinh(b))^2+(cosh(a)^2*sinh(b)+2*sinh(a)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b\
))+1))*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(a)*%e^(lambert_w(-%e^-1*sinh(a)*\
sinh(b))+1)+2*cosh(a)^2*sinh(b))/(cosh(a)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b)\
)+1)+2*cosh(a)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(\
a)*sinh(b))+cosh(a)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*\
sinh(a)*sinh(b))^2) 
0 2 2 -(cosh(a)*sinh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
0 3 3 -(cosh(a)*sinh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
0 4 4 -(cosh(a)*sinh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
1 1 1 (sinh(b)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(\
a)*sinh(b))^2+(sinh(a)*cosh(b)^2+2*sinh(b)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b\
))+1))*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(b)*%e^(lambert_w(-%e^-1*sinh(a)*\
sinh(b))+1)+2*sinh(a)*cosh(b)^2)/(cosh(b)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b)\
)+1)+2*cosh(b)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(\
a)*sinh(b))+cosh(b)*%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*\
sinh(a)*sinh(b))^2) 
1 2 2 -(sinh(a)*cosh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
1 3 3 -(sinh(a)*cosh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
1 4 4 -(sinh(a)*cosh(b))/(%e^(lambert_w(-%e^-1*sinh(a)*sinh(b))+1)+2*%e^(lambe\
rt_w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))+%e^(lambert_\
w(-%e^-1*sinh(a)*sinh(b))+1)*lambert_w(-%e^-1*sinh(a)*sinh(b))^2) 
2 2 0 -(sinh(a)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
2 2 1 -(sinh(b)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
3 3 0 -(sinh(a)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
3 3 1 -(sinh(b)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
4 4 0 -(sinh(a)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(a))/(2*cosh(a)) 
4 4 1 -(sinh(b)*lambert_w(-%e^-1*sinh(a)*sinh(b))+sinh(b))/(2*cosh(b)) 
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
    c.gsub!("%e^-1","exp(-1)")
    c.gsub!("%e^a","exp(a)")
    c.gsub!("%e^b","exp(b)")
    c.gsub!("%e^-a","exp(-a)")
    c.gsub!("%e^-b","exp(-b)")
    c.gsub!("%e^-lambert_w(-(sinh(a)*sinh(b))/e)","exp(-lambert_w(-(sinh(a)*sinh(b))/e))")
    c.gsub!("%e^","exp")
    c.gsub!("^","**")
    c.gsub!("%e","math.e")
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
    c.gsub!("%e^-1","exp(-1)")
    c.gsub!("%e^a","e^a ")
    c.gsub!("%e^b","a^b ")
    c.gsub!("%e^-a","e^{-a}")
    c.gsub!("%e^-b","e^{-b}")
    c.gsub!("%e^","exp")
    #c.gsub!("lambert_w(exp((-b)+a-1))","W(e^{-b+a-1})")
    c.gsub!("lambert_w(exp((-b)+a-1))","Q")
    c.gsub!("\*","")
    c.gsub!("lambert_w","W")
    c.gsub!("W(-exp(-1)sinh(a)sinh(b))+1","r")
    c.gsub!("W(-exp(-1)sinh(a)sinh(b))","(r-1)")
    c.gsub!("%e^","e^")
    c.gsub!(/(exp|cosh|sinh)/) {"\\#{$1}"}
    print "\\Gamma\\indices{^#{coord[k]}_#{coord[i]}_#{coord[j]}} &= #{c} \\\\\n"
  end
end

end
