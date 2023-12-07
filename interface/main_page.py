import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd
import math
from Notebook import *
#from running import *

class main_page(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')

        titre_text = "MAIN"

        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

        self.values = []

        self.texts = ["Model Param","PID propeller","PID controller","od Cycle","MVEM model","Param energy","Excell imports"]

        self.notebook = ttk.Notebook(self, bootstyle="info")

        self.frames_notebook = [Notebook_Model_param(self.notebook),
                                Notebook_propeller(self.notebook),
                                Notebook_PIDcontroller(self.notebook),
                                Notebook_0dCycle(self.notebook),Notebook_MVEM_model(self.notebook),
                                Notebook_Param_energy(self.notebook),Notebook_excell_sheets(self.notebook)]
        
        for count,text in enumerate(self.texts):
            self.notebook.add(self.frames_notebook[count],text=text)

    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.notebook.pack(fill=BOTH,expand=YES)
        self.update_idletasks()


    def unshow(self):
        self.pack_forget()

    def return_all_values(self):

        self.values = []

        for i in range(len(self.frames_notebook)-1):
            self.values.extend(self.frames_notebook[i].return_values())

        return self.values
    
    def return_file_paths(self):

        return self.frames_notebook[6].return_values()
    
   # def run_page(self):
    #   Run_program()
    


        
        

