#!/usr/bin/ruby

# For documentation, see the README.md file in the same subdirectory as this script.

require 'digest'
require 'fileutils'

def main()
  t = gets(nil)
  if t.nil? then exit(-1) end
  $saved_comments = Hash.new
  $vars_placeholders = Hash.new
  t = process(t)
  print t
end

def process(t)
  t = preprocess(t)
  t = translate_stuff(t)
  t = postprocess(t)
  return t
end

def translate_stuff(t)
  t.gsub!(/^(\s*)def\s+([^:]+):/) {$1+"function "+$2}
  t.gsub!(/^(\s*)if\s+([^:]+):/) {$1+"if ("+$2+")"}
  t.gsub!(/^(\s*)else:/) {$1+"else"}
  t.gsub!(/^(\s*)(from|import).*/) {''}
  t.gsub!(/ and /,' && ')
  t.gsub!(/ or /,' || ')
  t.gsub!(/ not /,' ! ')
  return t
end

def postprocess(t)
  $saved_comments.each_pair { |key,value|
    t.gsub!(/#{key}/,value)
  }
  $vars_placeholders.each_pair { |key,value|
    vars = value.keys
    if vars.length>0 then
      decl = "var "+vars.join(',')+";\n"
    else
      decl = 'wugga'
    end
    t.gsub!(/#{key}/,decl)
  }
  t.gsub!(/ \/\*\s*"""/,"/*") # kludge for docstrings
  t.gsub!(/^((  )*) \/\*/) {$1+"/*"} # kludge to convert 2n+1 spaces before /* to 2n
  return t
end

def preprocess(t)
  t = docstrings_to_comments(t)
  # Protect text of comments from munging:
  t.gsub!(/#([^\n]*)$/) {d=digest($1); $saved_comments[d]=$1; '#'+digest($1); }
  # Combine long lines into single lines:
  t.gsub!(/([^\n]*)\\\s*\n([^\n]*)/) {"#{$1} #{$2}"}
  # Change tabs to 8 blanks:
  t.gsub!(/\t/,"        ")
  # Hand-translated lines are marked with #js comments.
  t.gsub!(/^( *)([^#\n]*)#js( *)([^\n]*)$/) {"#{$1}#{$4}"}
  # Translate comments:
  t.gsub!(/^([^#]*)#([^\n]*)$/) {"#{$1}/*#{$2} */"}
  # Split into lines:
  lines = t.split(/\n+/)
  # Curly braces around indented blocks:
  t2 = ''
  lines.push("___DELETEME___") # to trigger }'s at end
  indent_stack = [0]  # level of indentation at the start of the block
  opener_stack = [''] # keyword that started the block, e.g. "def"
  last_line = ''
  current_function_key = nil
  0.upto(lines.length-1) { |i|
    l = lines[i]
    l =~ /^( *)/
    ind = $1.length
    done = false
    if ind>indent_stack[-1] then
      t2.gsub!(/;\n\Z/,'') # go back and erase semicolon and newline from previous line
      t2 = t2 + " {\n" + l + ";\n"
      indent_stack.push(ind)
      last_line =~ /([^ ]*)/ # first non-blank text is assumed to be keyword starting the block
      opener = $1
      opener_stack.push(opener)
      if opener=='def' then
        current_function_key = "vars_placeholder"+digest(l+i.to_s)
        t2 = t2+' '*ind+current_function_key
        $vars_placeholders[current_function_key] = Hash.new
      end
      done = true
    end
    while ind<indent_stack[-1]
      indent_stack.pop
      opener = opener_stack.pop
      if opener=='def' || opener=='' then semicolon='' else semicolon=';' end
      t2 = t2 + (" "*indent_stack[-1]) + "}" + semicolon +"\n"
    end
    if !done then
      if l=~/^\s*(\w[\w0-9,\[\]]*)\s*=/ then l=translate_assignment(l,current_function_key) end
      t2 = t2 + l + ";\n"
    end
    last_line = l
  }
  t = t2
  t.sub!(/___DELETEME___;?\n?/,'')
  # Bring semicolons and curly braces to the left of comments:
  t.gsub!(/(\/\*[^\n]*\*\/);$/) {"; #{$1}"}
  t.gsub!(/(\/\*[^\n]*\*\/) *{$/) {"{ #{$1}"}
  # No semicolon on a line consisting only of a comment:
  t.gsub!(/^( *);( *)(\/\*[^\n]*\*\/)( *)$/) {"#{$1+$2+$3+$4}"}
  # Done:
  return t
end

def translate_assignment(l,current_function_key)
  l=~/^\s*(.*)\s*=(.*)/
  lhs,rhs = $1,$2
  if lhs.nil? || rhs.nil? then return l end
  lhs.scan(/[a-zA-Z]\w*/) { |var| 
    # In multiple assignment like x,y,z=array, loop over variables. In something like g[i]=...,
    # this also has the effect of declaring i, which is harmless.
    if !current_function_key.nil? then
      $vars_placeholders[current_function_key][var] = 1
    end
  }
  has_math = (rhs=~/(sin|cos|tan|exp|log|sqrt|abs|\*\*)/)
  if has_math then
    rhs.gsub!(/\*\*/,'^') # translate exponentiation to maxima syntax
    rhs.gsub!(/^\s+/,'') # delete leading whitespace, which upsets it, not sure why
    if rhs=~/(.*[^\s])\s*\/\*(.*)\*\// then
      rhs,comment = $1,$2
    else
      comment = ''
    end
    file1 = "temp1_"+digest(rhs)+".mac"
    file2 = "temp2_"+digest(rhs)+".mac"
    File.open(file1,'w') { |f|
      f.print rhs
    }
    cmd = "python3 translate_maxima.py <#{file1} >#{file2}"
    err = ''
    if !system(cmd) then
      err = "cmd=#{cmd} failed, $?=#{$?}\n"
    end
    if err=='' then
      File.open(file2,'r') { |f|
        result = f.gets(nil)
        if ! result.nil? then 
          if result=~/^\^(.*)/ then
            err = $1
          else
            l=lhs+"="+result
          end
        else
          err = "null result"
        end
      }
    end
    remove_temp_file(file1)
    remove_temp_file(file2)
    if err!='' then comment=err+' '+comment end
    if comment!='' then l=l+' /*'+comment+'*/' end
  end
  return l
end

def remove_temp_file(file)
  if FileTest.exist?(file) then FileUtils.rm(file) end
end

def docstrings_to_comments(t)
  lines = t.split(/\n+/)
  inside_docstring = false
  t2 = ''
  inside = false
  k = 1
  this_comment = ''
  lines.each { |l|
    done = false
    if l=~/( *)"""(.*)/ then
      done = true
      if not inside then
        this_comment = $2
        t2 = t2+$1+"/*"
        
      else
        this_comment = this_comment+$1
        d = digest(k.to_s)
        k = k+1
        $saved_comments[d]=this_comment
        this_comment = ''
        l = d+"*/"+$2
      end
      inside = !inside
    end
    if true then
      if !inside then
        t2 = t2+l+"\n"
      else
        this_comment = this_comment+l+"\n"
      end
    end
  }
  t2 = t2+"\n"
  return t2
end

def digest(s)
  return Digest::SHA256.hexdigest(s).upcase
  # upcase makes it less likely that something inside the hex will match a keyword or something
end

main()
