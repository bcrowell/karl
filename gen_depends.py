#!/usr/bin/python3

# This script generates the file "depend," which is included in the Makefile.
# To run this script, do a "make depend."

import re

def main():
  targets = {}

  generic = {"test":False}

  physics_and_math = {**generic,**{"name":"physics_and_math"}}
  for target in prep_list("""
        math_util      
        lambert_w_stuff
        schwarzschild  
        kruskal        
        keplerian
        transform      
        runge_kutta
        fancy
        angular        
        vector
        conserved
        euclidean
        celestial
      """):
    add_target(targets,target,physics_and_math)

  util = {**generic,**{"name":"util"}}
  for target in prep_list("""
        io_util        
        c_libs
        test
      """):
    add_target(targets,target,util)

  test = {**generic,**{"name":"test","test":True}}
  for target in prep_list("""
        lambert_w 
        schwarzschild  
        kruskal
        keplerian
        transform
        runge_kutta
        fancy
        angular
        math_util
        math
        euclidean
        celestial
      """):
    add_target(targets,target,test)

  print("# --------------------------------------------------------------------------------")
  print("# This file was generated automatically by gen_depends.py. Do not edit it by hand.")
  print("# --------------------------------------------------------------------------------")
  all_tests = ""
  all_test_obj = []
  all_js = []
  all_py = []
  all_js_modules = []
  all_py_modules = []
  individual = ""
  include_files = "src/include/*.h"
  for npass in range(2):
    if npass==1:
      js_modules = ' '.join(all_js_modules)
      py_modules = ' '.join(all_py_modules)
      print("test: "+' '.join(all_test_obj)+' '+js_modules+' '+py_modules+"\n"+all_tests)
      print("js: "+' '.join(all_js)+"\n"+"	#")
      print("py: "+' '.join(all_py)+"\n"+"	#")
    for pat_name in targets:
      for target in targets[pat_name]:
        pat = targets[pat_name][target]
        is_test = ("test" in pat and pat["test"])
        if is_test:
          prefix = "test_"
          src_dir = "test/"
        else:
          prefix = ""
          src_dir = "physics/"
        pp = "src/"+src_dir+prefix+target+".pp"
        py = "obj/"+prefix+target+".py"
        js = "js/"+prefix+target+".js"
        jsi = "js/"+prefix+target+".jsi"
        if npass==0:
          all_js.append(js)
          all_py.append(py)
          if not is_test:
            all_js_modules.append(js)
            all_py_modules.append(py)
          if is_test:
            all_tests = all_tests + "	@make test_"+target+"\n"
            all_test_obj.append(py)
            all_test_obj.append(js)
        if npass==1:
          print("#--------- "+target+" ("+pat_name+")")
          print(py+": "+pp+" "+include_files)
          print("	filepp $(FILEPP_OPTIONS) -DLANG=python "+pp+" -o "+py)
          print("	@chmod +x "+py)
          print(js+": "+pp+" "+include_files)
          print("	filepp $(FILEPP_OPTIONS) -DLANG=js "+pp+" -o "+jsi)
          print("	pj/pj.rb "+jsi+" karl <"+jsi+" >"+js)
          print("	@rm "+jsi)
          print("	@-js-beautify --replace -n -s 2 "+js)
          if is_test:
            print("test_"+target+"_py: "+py+" "+py_modules+" obj/karl.so")
            print("	@$(PYTHON3) "+py)
            print("test_"+target+"_js: "+js+" "+js_modules)
            print("	@$(PYTHON3) run_js.py "+js)
            print("test_"+target+":")
            print("	@make -s --no-print-directory test_"+target+"_py")
            print("	@make -s --no-print-directory test_"+target+"_js")
          print("")

def add_target(targets,target,pat):
  pat_name = pat["name"]
  if not (pat_name in targets):
    targets[pat_name] = {}
  targets[pat_name][target] = pat

def prep_list(s):
  a = re.split("\s+",s)
  a = map(lambda x : re.sub(r"\s",'',x), a)
  a = map(lambda x : re.sub(r"#.*",'',x), a)
  a = filter(lambda x : x!='', a)
  return a

main()
