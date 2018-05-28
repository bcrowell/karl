#!/usr/bin/ruby

# This is to clean up maxima's output of the christoffel coefficients for K-S coordinates.
# There is output for python and latex. Can comment out one or the other as desired.
# The symbol om corresponds to what I'm now calling \ell, or r-1.
# Order of indices is ctensor's:
#    symmetric on 1st 2 indices
#    contravariant on final index

# (%o3)                         [v, w, theta, phi]

stuff = <<STUFF;
1 1 1 (-w/(om*%e^(om+1)+2*om^2*%e^(om+1)+om^3*%e^(om+1)))-om^2/(v+2*om*v+om^2*v)-(2*om)/(v+2*om*v+om^2*v)-1/(v+2*om*v+om^2*v) 
1 3 3 -w/(%e^(om+1)+2*om*%e^(om+1)+om^2*%e^(om+1)) 
1 4 4 -w/(%e^(om+1)+2*om*%e^(om+1)+om^2*%e^(om+1)) 
2 2 2 (-om^2/(w+2*om*w+om^2*w))-(2*om)/(w+2*om*w+om^2*w)-1/(w+2*om*w+om^2*w)-v/(om*%e^(om+1)+2*om^2*%e^(om+1)+om^3*%e^(om+1)) 
2 3 3 -v/(%e^(om+1)+2*om*%e^(om+1)+om^2*%e^(om+1)) 
2 4 4 -v/(%e^(om+1)+2*om*%e^(om+1)+om^2*%e^(om+1)) 
3 3 1 (%e^((-om)-1)*v^2*w)/(8*om)+(%e^((-om)-1)*v^2*w)/8 
3 3 2 (%e^((-om)-1)*v*w^2)/(8*om)+(%e^((-om)-1)*v*w^2)/8 
3 4 4 cos(theta)/sin(theta) 
4 4 1 (%e^((-om)-1)*sin(theta)^2*v^2*w)/(8*om)+(%e^((-om)-1)*sin(theta)^2*v^2*w)/8 
4 4 2 (%e^((-om)-1)*sin(theta)^2*v*w^2)/(8*om)+(%e^((-om)-1)*sin(theta)^2*v*w^2)/8 
4 4 3 -cos(theta)*sin(theta) 
STUFF


###################################################################
# python output
###################################################################
coord = ['v','w','theta','phi']

stuff.each_line do |line|
  if line=~/(\d) (\d) (\d) (.*)/ then
    i,j,k,c = $1.to_i,$2.to_i,$3.to_i,$4
    c.gsub!("om","___omega___")
    c.gsub!("%e^","exp")
    c.gsub!("^","**")
    c.gsub!("___omega___","(r-1)")
    if i!=j then top=2 else top=1 end
    1.upto(top) { |m|
      if m==2 then u=i; i=j; j=u end
      print "  ch[#{i-1}][#{j-1}][#{k-1}] = #{c}\n  #   ... ^#{coord[k-1]} _#{coord[i-1]} #{coord[j-1]}\n"
    }
  end
end


###################################################################
# latex output
###################################################################

if false then

coord = ['v','w','\\theta','\\phi']

stuff.each_line do |line|
  if line=~/(\d) (\d) (\d) (.*)/ then
    i,j,k,c = $1.to_i,$2.to_i,$3.to_i,$4
    c.gsub!("om","___omega___")
    c.gsub!("sin(theta)^2","sin^2theta")
    c.gsub!("sin(theta)","sintheta")
    c.gsub!("cos(theta)","costheta")
    c.gsub!("sin","\\sin")
    c.gsub!("cos","\\cos")
    c.gsub!("theta","\\theta ")
    c.gsub!("phi","\\phi")
    c.gsub!("(-___omega___)-1","-___omega___-1")
    c.gsub!("%e^(-___omega___-1)","e^{-r}")
    c.gsub!("%e^(___omega___+1)","e^r")
    c.gsub!("___omega___","\\omega ")
    c.gsub!("\*","")
    print "\\Gamma\\indices{^#{coord[k-1]}_#{coord[i-1]}_#{coord[j-1]}} &= #{c} \\\\\n"
  end
end

end
