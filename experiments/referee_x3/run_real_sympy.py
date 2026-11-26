"""run a bundle script with real sympy in place of sympy_pari_shim (the one-line change, made by a wrapper)"""
import sys, os, runpy, sympy
sys.modules['sympy_pari_shim'] = sympy
script = os.path.abspath(sys.argv[1]); sys.argv = [script] + sys.argv[2:]
sys.path.insert(0, os.path.dirname(script))
runpy.run_path(script, run_name='__main__')
