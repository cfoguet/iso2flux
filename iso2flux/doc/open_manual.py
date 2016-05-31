import os
import sys

def open_manual():
   path = os.path.dirname(os.path.abspath(__file__))+"/manual.pdf"
   print path
   if sys.platform.startswith('win'):
      os.system("start "+path)
   elif sys.platform.startswith("linux"):
        os.system("xdg-open "+'"'+path+'"')
   else:
        os.system("open "+path)
