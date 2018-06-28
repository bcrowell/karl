#!/usr/bin/python3

import sys,os

def main():
  js = sys.argv[1]
  os.system("cd js ; rhino -opt -1 ../"+js+" ; cd -")

main()
