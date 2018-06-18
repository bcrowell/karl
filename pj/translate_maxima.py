#!/usr/bin/python3

# For documentation, see the README.md file in the same subdirectory as this script.
# Read lines containing expressions as maxima source code from stdin and translate them one at a time.

import os,os.path,sys,hashlib,re,copy

##################################################################################################
# The following global variables contain all the information about how to translate from Maxima to
# the target language:

operators = {'MTIMES':'*','MPLUS':'+','RAT':'/'}
# Operators in Maxima. The associative ones aren't just 2-valent, can have 3 or more inputs.
# There are no subtraction and division operators, they're handled as a+(-1*b) and a*b^-1.
# Although MEXPT (exponentiation) is an operator, if you want it to be changed into a function like
# pow(), just leave it off the list, and put it into fns_to_translate.

constants = {'%e':'math.E','%pi':'math.PI'}
fns_to_prepend_with_math = ['sqrt','abs','sin','cos','tan','sinh','cosh','tanh',
                            'arcsinh','arccosh','arctanh','exp','log',
                            'pow']
fns_to_translate = {'mexpt':'pow','mabs':'abs'}
##################################################################################################

def main():
  for line in sys.stdin:
    line = line.rstrip('\n') # strip newline
    if not re.match("[^\s]",line): continue # skip blank lines
    lisp = maxima_to_lisp(line) # will not return if there's an error
    lisp_lines = lisp.split("\n")
    # Remove first, second, and last lines, which are of no interest:
    lisp_lines.pop(0)
    lisp_lines.pop(0)
    lisp_lines.pop
    # Put the rest back together without newlines:
    lisp = " ".join(lisp_lines)
    lisp = re.sub(r"\|",'',lisp) # causes readlisp to go into an infinite loop: https://github.com/olemb/readlisp/issues/2
    #print("lisp=",lisp)
    print(translate(readlisp(lisp)))

def maxima_to_lisp(maxima_code):
  debug = False
  d = hashlib.md5(maxima_code.encode('utf-8')).hexdigest()
  temp1 = "/tmp/translate_maxima_"+d+".lisp"
  temp2 = "/tmp/translate_maxima_2_"+d
  temp3 = "/tmp/translate_maxima_3_"+d
  q = '\\"'
  # ... unix only, will typically be deleted on reboot even if our cleanup attempt below fails
  code = "myexpression:"+maxima_code+"; save("+q+temp1+q+",myexpression); quit();"
  cmd = "maxima -q --batch-string=\""+code+"\" 1>"+temp2+" 2>"+temp3
  if debug:
    print("temp1=",temp1,", code=",code,", cmd=",cmd)
  try:
    os.system(cmd)
    if os.path.isfile(temp1):
      with open(temp1, 'r') as lisp_file:
        lisp = lisp_file.read()
      if debug:
        print("lisp=",lisp)
    else:
      if debug:
        print('temp file not created?')
        os.system("cat "+temp2)
        os.system("cat "+temp3)
      exit(-1)
  finally:
    remove_file_if_it_exists(temp1)
  return lisp

def remove_file_if_it_exists(filename):
  if os.path.isfile(filename):
    os.remove(filename)

########################################################################################
# Translate:
########################################################################################
def translate(e):
  #print(e)
  e = clean_up(copy.deepcopy(e))
  try:
    return decompile(copy.deepcopy(e))
  except:
    return "^error processing "+str(e)

# Typical input:
#     ['MTIMES', ['MPLUS', 'a', 'b'], 'c']
def decompile(e):
  e = copy.deepcopy(e)
  if isinstance(e,list):
    if e[0] in operators:
      operator = e[0]
      e.pop(0)
      stuff = operators[operator].join(["("+decompile(x)+")" for x in e])
      return "("+stuff+")"
    if re.match("^%?[a-zA-Z_]+",e[0]):
      # A function, such as %SIN, or possibly some function we don't know about, wahoo(x) (no % sign).
      # May have more than one input
      f = re.sub("^%",'',e[0]).lower()
      if f in fns_to_translate:
        f = fns_to_translate[f]
      if f in fns_to_prepend_with_math:
        f = "math."+f
      e.pop(0)
      stuff = ",".join(["("+decompile(x)+")" for x in e])
      return f+"("+stuff+")"
    # by default, just return a list:
    return " ".join([decompile(x) for x in e])
  else:
    if e in constants: # a constant such as %PI
      return constants[e]
    return str(e)

########################################################################################
# Clean up the lisp input:
########################################################################################

# Lisp expression for (a+b)*c:
#
#    (DSKSETQ $MYEXPRESSION '((MTIMES SIMP) ((MPLUS SIMP) $A $B) $C)) 
#
# As python object:
#     [LispSymbol('DSKSETQ'), LispSymbol('$MYEXPRESSION'), LispSymbol("'"),
#          [
#            [LispSymbol('MTIMES'), LispSymbol('SIMP')], 
#            [
#              [LispSymbol('MPLUS'), LispSymbol('SIMP')],
#              LispSymbol('$A'),
#              LispSymbol('$B')
#            ],
#          LispSymbol('$C')]
#     ]
# After I do symbols_to_strings, this becomes:
#     ['DSKSETQ', '$MYEXPRESSION', "'", [['MTIMES', 'SIMP'], [['MPLUS', 'SIMP'], '$A', '$B'], '$C']]
# The ' is s-expression syntactic sugar for the quote operator, ' foo -> (quote foo),
# basically prevents evaluation. SIMP means the expression has already been simplified.
# After further processing:
#     ['MTIMES', ['MPLUS', 'a', 'b'], 'c']

def clean_up(e):
  # Translate all the LispSymbol objects into strings, and remove the initial
  # remove 'DSKSETQ', '$MYEXPRESSION', "'" :
  debug = False
  if debug: print(e)
  e = copy.deepcopy(e)
  if not isinstance(e[2],list) and str(e[2])=="'": e.pop(2)
  e = copy.deepcopy(e)
  e = symbols_to_strings(copy.deepcopy(e))[2]
  if debug: print(e)
  e = remove_simp(copy.deepcopy(e))
  if debug: print(e)
  e = remove_uppercase_and_dollar_signs(copy.deepcopy(e))
  if debug: print(e)
  return e

def symbols_to_strings(e):
  if isinstance(e,LispSymbol):
    return str(e)
  if isinstance(e,float):
    return str(e)
  if isinstance(e,list):
    return [symbols_to_strings(x) for x in e]

def remove_simp(e):
  if isinstance(e,str): 
    return e
  if isinstance(e,list) and len(e)==2 and e[1]=='SIMP':
    return e[0]
  else:
    return [remove_simp(x) for x in e]

def remove_uppercase_and_dollar_signs(e):
  if isinstance(e,str):
    e = re.sub(r"^'",'',e) # remove ' from stuff like '$%pi
    if re.match(r"\$",e): # variable name
      e = e.lower() # lowercase
      e = re.sub(r"^\$",'',e) # remove $ from symbols
    return e
  if isinstance(e,list):
    return [remove_uppercase_and_dollar_signs(x) for x in e]

########################################################################################
# Ole Martin Bjorndalen's lisp parser
# https://github.com/olemb/readlisp
########################################################################################
"""
readlisp 0.2 - a lisp parser

Ole Martin Bjorndalen
ombdalen@gmail.com
http://nerdly.info/ole/

License: MIT

2006-12-21
"""

import io
import types

EOF = ''
END_PAREN = ')'

class CharFile:
    def __init__(self, file):
        self.file = file
        self.c = ''

    def getchar(self):
        if self.c:
            c = self.c
            self.c = ''
            return c
        else:
            self.c = ''
            c = self.file.read(1)
            # print [c]
            return c

    def ungetchar(self, c):
        self.c = c

class LispSymbol:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return 'LispSymbol(%s)' % repr(self.name)

    def __str__(self):
        return self.name

class LispReader:
    def __init__(self, file):
        self.file = CharFile(file)
        self.whitespace = ' \r\n\t'
        self.atom_end = self.whitespace + '()"|;'

    def _skip_whitespace(self):
        """Read until next non-whitespace or EOF"""
        while 1:
            c = self.file.getchar()
            if c == EOF or c not in self.whitespace:
                self.file.ungetchar(c)
                return

    def _skip_comment(self):
        """Read until EOL or EOF"""
        while 1:
            c = self.file.getchar()
            if c == '\n' or c == EOF:
                return

    def _parse_atom(self, atom):
        """Parse an atom and return and integer, a float or a symbol"""
        try:
            return int(atom)
        except ValueError:
            try:
                return float(atom)
            except ValueError:
                return LispSymbol(atom)

    def _read_atom(self):
        """Read a symbol or number"""
        atom = ''

        while 1:
            c = self.file.getchar()
            if c in self.atom_end or c == EOF:
                self.file.ungetchar(c)
                break
            atom += c

        return self._parse_atom(atom)

    def _read_string(self):
        string = ''
        
        while 1:
            c = self.file.getchar()

            if c == '':
                raise EOFError('Missing end quote for string')
            elif c == '"':
                return string
            elif c == '\\':
                c = self.file.getchar()
                if c == '':
                    raise EOFError('EOF after escape sequence')
                else:
                    string += c
            else:
                string += c

    def _read_list(self):        
        items = []
        while 1:
            item = self._read_expr()
            if item is END_PAREN:
                break
            else:
                items.append(item)
        return items

    def _read_expr(self):

        while 1:
            self._skip_whitespace()

            c = self.file.getchar()

            if c == EOF:
                raise EOFError('End of LISP stream')
            elif c == ';':
                self._skip_comment()
            elif c == '(':
                return self._read_list()
            elif c == ')':
                return END_PAREN
            elif c == '"':
                return self._read_string()
            #elif c == '|':
            #    return self.read_quoted_symbol()
            else:
                self.file.ungetchar(c)
                return self._read_atom()

    def __iter__(self):
        while 1:
            try:
                yield self._read_expr()
            except EOFError:
                return


def writelisp(obj):
    """Convert a python object into an equivalent lisp expression."""

    if type(obj) is types.ListType:
        return '(%s)' % ' '.join(map(writelisp, obj))
    elif type(obj) is types.StringType:
        out = '"'
        for c in obj:
            if c in '\\"':
                out += '\\'
            out += c
        out += '"'
        return out
    elif type(obj) in [types.LongType, types.IntType]:
        return str(obj)
    elif type(obj) is types.ComplexType:
        return '#C(%s %s)' % (obj.real, obj.imag)
    elif obj == None:
        return 'nil'
    else:
        return repr(obj)

def readlisp(text):
    """Read the first lisp expression in the string"""
    return LispReader(io.StringIO(text))._read_expr() # modified for python3 compatibility - 2018 jun 17, BC

########################################################################################
# main
########################################################################################

if __name__ == '__main__':
  main()
