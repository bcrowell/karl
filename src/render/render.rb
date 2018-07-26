#!/usr/bin/ruby

require 'csv'
require 'chunky_png'
  # ubuntu package ruby-chunky-png
require 'hsluv' 
  # http://www.hsluv.org
  # https://github.com/hsluv/hsluv-ruby
  # sudo gem install hsluv

$PI=3.1415926535

$gamma = 0.75 # https://en.wikipedia.org/wiki/Gamma_correction
# A value less than 1 makes stars appear more uniform in brightness, makes
# dim stars easier to see.

$default_width = 1000

def main
  blur = 1 # std dev. of gaussian blur, in units of pixels
  truncate_blur = 4.0 # # truncate blur at 3 s.d.
  overexpose = 10.0
  w = $default_width
  image_i = Array.new(w) { Array.new(w) } # intensity
  image_h = Array.new(w) { Array.new(w) } # intensity-averaged hue
  image_s = Array.new(w) { Array.new(w) } # intensity-averaged saturation
  0.upto(w-1) { |i|
    0.upto(w-1) { |j|
      image_i[i][j] = 0.0
      image_h[i][j] = 0.0
      image_s[i][j] = 0.0
    }
  }
  CSV.foreach("stars.csv") do |row|
    alpha,phi,brightness,bv = row[0].to_f,row[1].to_f,row[2].to_f,row[3].to_f
    if false then
      r = w*alpha/($PI/2.0)
      # ... rough and ready projection; away from b.h.; gives roughly 45-degree radius field of view
    else
      r = w*($PI-alpha)/($PI)
      # towards b.h., bigger field of view
    end
    x = 0.5*w+r*Math::cos(phi)
    y = 0.5*w-r*Math::sin(phi)
    next unless x>=0 && x<=w-1 && y>=0 && y<=w-1
    p = (blur*truncate_blur+0.5).to_i # size of square and inscribed circle on which to do computations
    blur2 = blur*blur
    (-p).upto(p) { |i|
      xx = x+i
      next unless xx>=0 && xx<=w-1
      (-p).upto(p) { |j|
        yy = y+j
        next unless yy>=0 && yy<=w-1
        # gaussian blur
        h = i*i+j*j
        next unless h/blur2<truncate_blur*truncate_blur # cut off to a circle, so it doesn't look boxy
        b = brightness*Math::exp(-0.5*h/blur2)
        image_i[xx][yy] += b
        hue,sat = bv_to_color(bv)
        image_h[xx][yy] += b*hue
        image_s[xx][yy] += b*sat
      }
    }
  end
  max = 0.0
  0.upto(w-1) { |i|
    0.upto(w-1) { |j|
      if image_i[i][j]>max then max=image_i[i][j] end
    }
  }
  if max==0.0 then $stderr.print("max=0"); exit(-1) end
  black = hsv_to_color(0,0,0)
  white = hsv_to_color(0,0,100.0)
  png_image = ChunkyPNG::Image.new(w,w,black)
  0.upto(w-1) { |i|
    0.upto(w-1) { |j|
      z = image_i[i][j]/max
      next if z==0.0
      z = z**$gamma
      l = 100.0*z
      l = l*overexpose # allow the brightest stars to overexpose the pixel
      if l>100.0 then l=100.0 end
      hue = put_in_range(image_h[i][j]/image_i[i][j],0.0,360.0)
      sat = put_in_range(image_s[i][j]/image_i[i][j],0.0,100.0)
      next if hue.nan? or sat.nan? or l.nan?
      png_image[i,j] = hsv_to_color(hue,sat,l)
    }
  }
  png_image.save("stars.png")
end

def put_in_range(x,min,max)
  if x<min then x=min end
  if x>max then x=max end
  return x
end

def bv_to_color(bv)
  # https://en.wikipedia.org/wiki/Color_index
  # The following is just a rough-and-ready approximation that I made up.
  # Better: http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/
  h = -((0.0-240.0)/(1.40-(-0.33)))*(bv-1.40)
  if h<0.0 then h=h+360.0 end
  if h>360.0 then h=h-360.0 end
  s = (bv-0.5)**4*100.0
  if s<0.0 then s=0.0 end
  if s>100.0 then s=100.0 end
  return [h,s]
end

def hsv_to_color(h,s,v)
  r,g,b = Hsluv::hsluv_to_rgb(h,s,v)
  return ChunkyPNG::Color::rgb((r*255).to_i,(g*255).to_i,(b*255).to_i)
end

main()
