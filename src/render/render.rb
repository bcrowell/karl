#!/usr/bin/ruby

require 'csv'
require 'chunky_png'
  # ubuntu package ruby-chunky-png
require 'hsluv' 
  # http://www.hsluv.org
  # https://github.com/hsluv/hsluv-ruby
  # sudo gem install hsluv
require 'json'

def main
  infile = 'stars.json'
  fov_deg = 45.0 # field of view in degrees
  view_rot_deg = 0.0 # angle away from black hole, i.e., 180 means looking directly away from it
  #----
  begin
    image_arrays = JSON.parse File.read(infile)
  rescue JSON::ParserError
    $stderr.print "invalid JSON syntax"
    exit(-1)
  end 
  #----
  image_i,image_h,image_s = image_arrays
  w = image_i.length # width of square field of view, in pixels
  #----
  black = hsv_to_color(0,0,0)
  white = hsv_to_color(0,0,100.0)
  png_image = ChunkyPNG::Image.new(w,w,black)
  max_hue = 265.0
  0.upto(w-1) { |i|
    0.upto(w-1) { |j|
      l = 100.0*image_i[i][j].to_f
      l = put_in_range(l,0.0,100.0)
      hue = put_in_range(image_h[i][j].to_f*max_hue,0.0,max_hue)
      sat = put_in_range(image_s[i][j].to_f*100.0,0.0,100.0)
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

def hsv_to_color(h,s,v)
  r,g,b = Hsluv::hsluv_to_rgb(h,s,v)
  return ChunkyPNG::Color::rgb((r*255).to_i,(g*255).to_i,(b*255).to_i)
end

main()
