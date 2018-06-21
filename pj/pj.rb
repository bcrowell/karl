#!/usr/bin/ruby

# For documentation, see the README.md file in the same subdirectory as this script.

require 'digest'
require 'fileutils'

$modules_not_to_load = ['numpy','scipy','math','sys']
# If the python code does an "import numpy", don't try to translate that into javascript.

def main()
  module_name = ''
  if ARGV.length>0 then module_name = extract_file_stem(ARGV[0]) end
  t = gets(nil)
  if t=~/^((\s*)if\s+(.*): *[^\s]+)/ then
    $stderr.print "Error: if statement on a single line: \'#{$1}\'\n"
    exit(-1)
  end
  if t.nil? then exit(-1) end
  is_main = false
  if t=~/\A#!/ then is_main=true end # If it has #!/usr/bin/python3 at the top, it's a program, not a module.
  $protected_strings = Hash.new
  $vars_placeholders = {"vars_placeholder_global"=>{}}
  t = process(t,module_name,is_main)
  print t
end

def extract_file_stem(filename)
  return filename.sub(/\A(.*)\//,'').sub(/\.[a-z]+/,'')
end

def process(t,module_name,is_main)
  t = preprocess(t)
  t = translate_stuff(t,module_name,is_main)
  t = "vars_placeholder_global\n"+t
  t = postprocess(t)
  t = header(module_name,is_main) + t
  return t
end

def header(module_name,is_main)
  if is_main then
    h = <<-HEADER
      /*
         --- main program #{module_name} ---
         This was translated from python. Do not edit directly.
      */
      karl = {};
      karl.modules_loaded = {};
      if (!(typeof window !== 'undefined')) { /* IS_BROWSER can't be defined yet. */
        load("lib/loader.js");
      }
      HEADER
  else
    h = <<-HEADER
      /*
         --- module #{module_name} ---
         This was translated from python. Do not edit directly.
      */
      if (typeof #{module_name} === 'undefined') {
        #{module_name} = {};
      }
      HEADER
  end
  return h
end

def translate_stuff(t,module_name,is_main)
  lines = t.split(/\n+/)
  t2 = ''
  0.upto(lines.length-1) { |i|
    l = lines[i]
    t2 = t2 + translate_stuff_one_line(l,module_name,is_main) + "\n"
  }
  return t2
end

def translate_stuff_one_line(t,module_name,is_main)
  # There is something wrong with the following, which are based on curlycurlycurly in footex. They
  # cause infinite loops, don't behave as expected. Don't use them without figuring out what's wrong.
  # paren = "(?:(?:\([^()]*\)|[^()]*)*)"; # match anything, as long as parens in it are matched, and not nested
  # parenparen = "(?:(?:\(#{paren}\)|#{paren})*)"; # allow one level of nesting
  # parenparenparen = "(?:(?:\(#{parenparen}\)|#{parenparen})*)"; # allow two levels of nesting

  # kludge: def, if, ... are handled here, but assignments and for loops are handled in preprocess()
  t.gsub!(/^(\s*)def\s+([^\(]+)([^:]+):/) {
    indentation,func,args = [$1,$2,$3];
    if is_main then m="" else m="#{module_name}." end
    "#{indentation}#{m}#{func} = function#{args}"
  }
  t.gsub!(/^(\s*)if\s+(.*):/) {$1+"if ("+$2+")"}
  t.gsub!(/^(\s*)else:/) {$1+"else"}
  t.gsub!(/^(\s*)import\s+([^;]*)/) {
    indentation,modules = [$1,$2]
    indentation+modules.split(/,/)\
      .select{ |m| !$modules_not_to_load.include?(m)}\
      .map{ |m| "karl.load(\"#{m}\")"}\
      .join(";")
  }
  t.gsub!(/^(\s*)from.*/) {''}
  t.gsub!(/(?<!\w)and(?!\w)/,' && ')
  t.gsub!(/(?<!\w)or(?!\w)/,' || ')
  t.gsub!(/(?<!\w)not(?!\w)/,' ! ')
  return t
end

def postprocess(t)
  t.gsub!(/__NO_SEMICOLON__;?/,'')
  $protected_strings.each_pair { |key,value|
    t.gsub!(/#{key}/,value)
  }
  $vars_placeholders.each_pair { |key,value|
    vars = value.keys
    if vars.length>0 then
      decl = "var "+vars.join(',')+";\n"
    else
      decl = ''
    end
    t.gsub!(/#{key}/,decl)
  }
  t.gsub!(/ \/\*\s*"""/,"/*") # kludge for docstrings
  t.gsub!(/^((  )*) \/\*/) {$1+"/*"} # kludge to convert 2n+1 spaces before /* to 2n
  t.gsub!(/__NO_TRANSLATION__/,'')
  return t
end

def preprocess(t)
  t = docstrings_to_comments(t)
  # Combine long lines into single lines:
  t.gsub!(/([^\n]*)\\\s*\n([^\n]*)/) {"#{$1} #{$2}"}
  # Change tabs to 8 blanks:
  t.gsub!(/\t/,"        ")
  # Hand-translated lines are marked with #js comments.
  t = hand_translated_lines(t)
  #$stderr.print "1"*80+"\n"+t # qwe
  # Protect text of comments from munging:
  t.gsub!(/#([^\n]*)$/) {d=digest("comment"+$1); $protected_strings[d]=$1; '#'+d }
  # Protect string literals:
  t.gsub!(/"([^"\n]*)"/) {d=digest("string_literal"+$1); $protected_strings[d]=$1; '"'+d+'"' }
  # Translate comments from # to /**/:
  t.gsub!(/^([^#]*)#([^\n]*)$/) {"#{$1}/*#{$2} */"}
  # Split into lines:
  lines = t.split(/\n+/)
  # Curly braces around indented blocks:
  t2 = ''
  lines.push("___DELETEME___") # to trigger }'s at end
  indent_stack = [0]  # level of indentation at the start of the block
  opener_stack = [''] # keyword that started the block, e.g. "def"
  last_line = ''
  current_function_key = "vars_placeholder_global"
  0.upto(lines.length-1) { |i|
    l = lines[i]
    l =~ /^( *)/
    ind = $1.length
    no_trans = (l=~/__NO_TRANSLATION__/)
    done = false
    if ind>indent_stack[-1] then
      done = true
      last_line =~ /([^ ]*)/ # first non-blank text is assumed to be keyword starting the block
      opener = process_opener($1)
      opener_stack.push(opener)
      if opener=='def' then
        current_function_key = "vars_placeholder"+digest(l+i.to_s)
      end
      indent_stack.push(ind)
      t2.gsub!(/;\n\Z/,'') # go back and erase semicolon and newline from previous line
      t2 = t2 + " {\n"
      if !no_trans then
        if opener=='def' then
          t2 = t2+' '*ind+current_function_key+"\n"
          $vars_placeholders[current_function_key] = Hash.new
        end
        t2 = t2 + translate_line(l,current_function_key) + ";\n"
      else
        t2 = t2+l
      end
    end
    while ind<indent_stack[-1]
      indent_stack.pop
      opener = opener_stack.pop
      t2 = t2 + (" "*indent_stack[-1]) + "}"
      if opener=='def' || opener=='' then t2=t2+"\n\n" else t2=t2+";\n" end
      if opener=='def' then current_function_key = "vars_placeholder_global" end
    end
    if !done then
      if !no_trans then
        t2 = t2 + translate_line(l,current_function_key) + ";\n"
      else
        t2 = t2 + l + "\n"
      end
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

def hand_translated_lines(t0)
  t = t0.clone
  lines = t.split(/\n+/)
  t2 = ''
  0.upto(lines.length-1) { |i|
    l = lines[i]
    t2 = t2 + hand_translated_line(l) + "\n"
  }
  return t2
end

def hand_translated_line(t0)
  t = t0.clone
  t.gsub!(/^( *)([^#\n]*)#js( *)([^\n]*)$/) {"#{$1}__NO_TRANSLATION__#{$4}__NO_SEMICOLON__"}
  return t
end

# Normally opener is whatever keyword opened the block, e.g., "def".
# But in a few cases, we get comments here. This happens, e.g., when include files have indented comments below
# defines.
def process_opener(opener)
  if opener=~/\/\*/ then return '' end
  return opener
end

# kludge: def, if, ... are handled in translate_stuff(), but assignments and for loops are handled in this
# routine, which gets called by preprocess()
def translate_line(l,current_function_key)
  if l=~/^\s*(\w[\w0-9,\[\]]*)\s*=/ then 
    return translate_assignment(l,current_function_key)
  end
  if l=~/for/ and l=~/range/ then
    return translate_for_loop(l,current_function_key)
  end
  if has_math(l) then l=translate_math_with_comment(l) end # things like function calls that have math in them
  return l
end

def translate_for_loop(l,current_function_key)
  # bug: doesn't preserve indentation, deletes comments; this should be handled using common logic, not cutting and
  #    pasting code from translate_assignment()
  # example: for j in range(3,len(x)):
  if l =~ /for\s+(\w+)\s+in\s+range\((.*)\):/ then
    var,r = [$1,$2]
    lo=0
    if r=~/(.*),(.*)/ then
      lo,hi = [$1,$2]
    else
      hi = r
    end
    return "for (var #{var}=#{lo}; #{var}<#{hi}; #{var}++)"
  end
end

# Heuristic pattern matching to see whether a line of code has math that needs to be
# translated.
def has_math(l)
  if (l=~/(sin|cos|tan|exp|log|sqrt|abs|lambert_w|floor|\*\*)/) then return true end
  return false
end

def translate_assignment(l,current_function_key)
  l=~/^(\s*)(.*)\s*=(.*)/
  indentation,lhs,rhs = [$1,$2,$3]
  if lhs.nil? || rhs.nil? then return l end
  list_variables_to_declare(lhs,current_function_key)
  if has_math(rhs) then
    rhs = translate_math_with_comment(rhs)
    l = indentation+lhs + "=" + rhs
  end
  l = translate_expression_common_idioms(l)
  return l
end

# This is very simpleminded, won't handle nested parens.
def translate_expression_common_idioms(e0)
  e = e0.clone
  # Translate len(x)->x.length, str(x)->x.toString().
  e.gsub!(/(?<!\w)len\(([^)]+)\)/) {"#{$1}.length"}
  e.gsub!(/(?<!\w)str\(([^)]+)\)/) {"#{$1}.toString()"}
  return e
end

# Input is like "y+z /* blah */", and the comment is maintained.
# Only handles the rhs of an assignment.
def translate_math_with_comment(e0)
  e = e0.clone
  if e=~/(.*[^\s])\s*\/\*(.*)\*\// then
    e,comment = $1,$2
  else
    comment = ''
  end
  e,err = translate_math_expression(e)
  if err!='' then comment=err+' '+comment end
  if comment!='' then e=e+' /*'+comment+'*/' end
  return e
end

# For a line like this
#     x = y+z[2] /* blah */,
# the input is "y+z[2]". If there is fancy math (basically just the exponentiation operator)
# that requires translate_maxima, use that. But in simpler cases, just do regexes, because
# (a) translate_maxima is super slow, and (b) translate_maxima reduces readability by adding
# lots of parens.
def translate_math_expression(e0)
  e = e0.clone
  err = ''
  if e=~/\*\*/ then # fancy math that requires translate_maxima
    e.gsub!(/\*\*/,'^') # translate exponentiation to maxima syntax
    e.gsub!(/^\s+/,'') # delete leading whitespace, which upsets it, not sure why
    # Protect array subscripts from translation, since translate_maxima can't handle them.
    # Also protect things like module.function(x), because maxima thinks . is multiplication.
    k = 0
    protect = Hash.new
    e.gsub!(/([a-zA-Z]\w*(\[[^\]]+\])+)/) {k=k+1; protect[k]=$1; "protectme#{k}"}
    e.gsub!(/(([a-zA-Z]\w*\.)+([a-zA-Z]\w*))/) {k=k+1; protect[k]=$1; "protectme#{k}"}
    e,err = translate_math_using_translate_maxima(e)
    e.gsub!(/protectme([0-9]+)/) {protect[$1.to_i]}
  else
    # Look for one of these math functions that is not preceded by a word char, and that is followed by a (.
    e.gsub!(/(?<!\w)(sqrt|abs|sin|cos|tan|sinh|cosh|tanh|arcsinh|arccosh|arctanh|exp|log|lambert_w|floor)(?=\()/) {
      "Math.#{$1}"
    }
    #     ... see similar list in fns_to_prepend_with_math in translate_maxima.
  end
  return [e,err]
end

def translate_math_using_translate_maxima(e0)
  e = e0.clone
  debug = false
  file1 = "temp1_"+digest(e)+".mac"
  file2 = "temp2_"+digest(e)+".mac"
  File.open(file1,'w') { |f|
    f.print e
  }
  my_dir = __dir__ # translate_maxima.py is supposed to be in my directory
  cmd = "python3 #{my_dir}/translate_maxima.py <#{file1} >#{file2}"
  err = ''
  if !system(cmd) then
    err = "cmd=#{cmd} failed, $?=#{$?}\n"
  end
  if debug then print "e=#{e}, ran command\n" end
  if err=='' then
    File.open(file2,'r') { |f|
      result = f.gets(nil).gsub(/\n$/,'')
      if debug then print "e=#{e} result=#{result}\n" end
      if ! result.nil? then 
        if result=~/^\^(.*)/ then
          err = $1
        else
          e = result
        end
      else
        err = "null result"
      end
    }
  end
  remove_temp_file(file1)
  remove_temp_file(file2)
  return [e,err]
end

def list_variables_to_declare(lhs,current_function_key)
  lhs.scan(/[a-zA-Z]\w*/) { |var|
    # Make a list of variables inside this function so that we can declare them.
    # In multiple assignment like x,y,z=array, loop over variables. In something like g[i]=...,
    # this also has the effect of declaring i, which is harmless.
    if !current_function_key.nil? then
      $vars_placeholders[current_function_key][var] = 1
    end
  }
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
        $protected_strings[d]=this_comment
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
  return Digest::SHA256.hexdigest(s).gsub(/[a-f]/,'')
  # Filter out alphabetic hex characters to prevent any parser from accidentally thinking this is something
  # other than an atomic symbol, not to be tampered with.
end

main()
