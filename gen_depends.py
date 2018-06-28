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
        transform      
        runge_kutta    
        angular        
        vector         
        test_math
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
        transform
        runge_kutta  
        angular
        math_util
        vector
      """):
    add_target(targets,target,test)

  print("# --------------------------------------------------------------------------------")
  print("# This file was generated automatically by gen_depends.py. Do not edit it by hand.")
  print("# --------------------------------------------------------------------------------")
  all_tests = ""
  all_test_obj = []
  all_js = []
  all_py = []
  for pat_name in targets:
    for target in targets[pat_name]:
      pat = targets[pat_name][target]
      is_test = ("test" in pat and pat["test"])
      prefix = ""
      if is_test:
        prefix = "test_"
      pp = "src/"+prefix+target+".pp"
      py = "obj/"+prefix+target+".py"
      js = "js/"+prefix+target+".js"
      jsi = "js/"+prefix+target+".jsi"
      include_files = "src/*.h"
      print("#--------- "+target+" ("+pat_name+")")
      print(py+": "+pp+" "+include_files)
      print("	filepp -DLANG=python "+pp+" -o "+py)
      print(js+": "+pp+" "+include_files)
      print("	filepp -DLANG=js "+pp+" -o "+jsi)
      print("	pj/pj.rb "+jsi+" karl <"+jsi+" >"+js)
      print("	@rm "+jsi)
      print("	@-js-beautify --replace -n -s 2 "+js)
      all_js.append(js)
      all_py.append(py)
      if is_test:
        print("test_"+target+": "+py+" "+js)
        print("	$(PYTHON3) "+py)
        print("	cd js ; rhino -opt -1 "+js+" ; cd -")
        all_tests = all_tests + "	make test_"+target+"\n"
        all_test_obj.append(py)
        all_test_obj.append(js)
      print()
  print("test: "+' '.join(all_test_obj)+"\n"+all_tests)
  print("js: "+' '.join(all_js)+"\n"+"	#")
  print("py: "+' '.join(all_py)+"\n"+"	#")

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
