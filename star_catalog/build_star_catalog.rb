#!/usr/bin/ruby

require 'sqlite3'
require "fileutils"

infile = '/usr/share/kstars/stars.dat'
outfile = 'mag7.sqlite'

def hms_to_radians(ss)
  if ss=~/(\d\d)(\d\d)(\d\d\.\d\d)/ then
    h,m,s = $1.to_f,$2.to_f,$3.to_f
    if h>24.0 || m>60 || s>60 then $stderr.print "error on input RA string #{ss}\n";exit(-1) end
    x = (h+m/60.0+s/3600.0)*Math::PI/12.0
    if x<0 || x>2.0*Math::PI then $stderr.print "error on input RA string #{ss}\n";exit(-1) end
    return x
  else
    $stderr.print "error on input RA string #{ss}\n"
    exit(-1)
  end
end

def dec_to_radians(ss)
  if ss=~/([\+\-]\d\d)(\d\d)(\d\d\.\d)/ then
    d,m,s = $1.to_f,$2.to_f,$3.to_f
    if d>90.0 || d<-90.0 || m>60 || s>60 then $stderr.print "error 1 on input dec string #{ss}\n";exit(-1) end
    x = (d+m/60.0+s/3600.0)*Math::PI/180.0
    if x<-0.5*Math::PI || x>0.5*Math::PI then $stderr.print "error 2 on input dec string #{ss}\n";exit(-1) end
    return x
  else
    $stderr.print "error 3 on input dec string #{ss}\n"
    exit(-1)
  end
end


FileUtils.rm_rf(outfile)
db = SQLite3::Database.new(outfile)
db.execute("CREATE TABLE stars(id INTEGER PRIMARY KEY, ra float, dec float, mag float, bv float)")
nlines = 0
File.readlines(infile).each { |line|
  if line =~ /^#/ then next end
  nlines = nlines+1
  if nlines%100==0 then print "." end
  if line=~/^(.{9,9}) (.{9,9}) [^ ]* (.{5,5})(.{5,5})/ then
    ra,dec,mag,bv = hms_to_radians($1),dec_to_radians($2),$3.to_f,$4.to_f
    if mag<7.0 then
      db.execute("INSERT INTO stars VALUES(#{nlines},#{ra},#{dec},#{mag},#{bv})")
    end
  end
}
print "\n"
db.close
$stderr.print "Successfully processed data records for #{nlines} stars.\n"


