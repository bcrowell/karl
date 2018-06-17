#!/usr/bin/ruby

# For documentation, see the README.md file in the same subdirectory as this script.

require 'digest'

def main()
  t = gets(nil)
  if t.nil? then exit(-1) end
  $saved_comments = Hash.new
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
  ['sin','cos','sqrt','abs','floor','log','exp'].each { |f|
    t.gsub!(/#{f}/) {"Math.#{f}"}
  }
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
  opener_stack = [''] # keyword that started the block, e.g. "function"
  last_line = ''
  0.upto(lines.length-1) { |i|
    l = lines[i]
    l =~ /^( *)/
    ind = $1.length
    done = false
    if ind>indent_stack[-1] then
      t2.gsub!(/;\n\Z/,'') # go back and erase semicolon and newline from previous line
      t2 = t2 + " {\n" + l + ";\n"
      indent_stack.push(ind)
      last_line =~ /([^ ]*)/ # first non-black text is assumed to be keyword starting the block
      opener_stack.push($1)
      done = true
    end
    while ind<indent_stack[-1]
      indent_stack.pop
      opener = opener_stack.pop
      if opener=='function' || opener=='' then semicolon='' else semicolon=';' end
      t2 = t2 + (" "*indent_stack[-1]) + "}" + semicolon +"\n"
    end
    if !done then
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
