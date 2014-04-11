import sys, compileall

sys.path.append(sys.argv[1])
compileall.compile_path()
