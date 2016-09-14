import Tkinter
import os
import tkFileDialog
import iso2flux
from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
from iso2flux.gui.gui import build_model_gui, launch_gui

def launch_build_model_gui():
    tk.destroy()
    global label_model
    root = Tkinter.Tk()
    gui = build_model_gui(root)
    root.mainloop()
    label_model=gui.label_model


def load_model():
    tk.withdraw()
    try:
      global label_model
      label_model=load_iso2flux_model(gui=True)
      tk.destroy()
    except:
       pass
    

if __name__ == "__main__": 
     tk=Tkinter.Tk()
     tk.withdraw()
     wd=tkFileDialog.askdirectory(title="Select working directory")
     os.chdir(wd)
     tk.deiconify()
     tk.title("Iso2flux")
     text=Tkinter.Label(tk,text="Welcome to Iso2Flux, select one option:   ")
     text.pack(side="top")
     load_button=Tkinter.Button(tk,text="Load Iso2Flux instance",command=load_model)
     create_button=Tkinter.Button(tk,text="Create Iso2Flux instance",command=launch_build_model_gui)
     load_button.pack()
     create_button.pack()
     tk.mainloop()
     launch_gui(label_model)
