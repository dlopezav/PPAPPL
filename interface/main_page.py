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

        self.frame = ttk.Frame(self, style='secondary')
        self.sauvegarde = ttk.Button(self.frame, command=lambda : self.sauvegarder_dans_excel(), style='success.Solid.TButton', text='sauvegarder',width=35) 
        self.recuperer = ttk.Button(self.frame, command=lambda : self.charger_depuis_excel(), style='danger.Solid.TButton', text='recuperer', width=35)
        self.executer = ttk.Button(self.frame, style='primary.Solid.TButton', text='executer', width=35, command=lambda: frameprincipal.show_values())
        

        self.frame.pack(side=BOTTOM)
        self.recuperer.pack(side=LEFT)
        self.sauvegarde.pack(side=LEFT)
        self.executer.pack(side=LEFT)

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
    
    def sauvegarder_dans_excel(self):
        nombre_archivo = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                                     filetypes=[("Excel Files", "*.xlsx")])
        # Création d'un classeur Excel
        if nombre_archivo:  # Si el usuario selecciona un archivo  
            wb = openpyxl.Workbook()
            sheet = wb.active
            sheet.title = "Datos"
            i = 1  # Empieza desde la primera fila
            for page in self.frames_notebook[:-1]:
                for j in range(len(page.float_vars)):
                    sheet.cell(row=i, column=2, value=page.float_vars[j].get())
                    i+=1
            wb.save(nombre_archivo)
            print(f"Données sauvegardées à : {nombre_archivo}")
    
    def charger_depuis_excel(self):
        nombre_archivo = filedialog.askopenfilename(filetypes=[("Excel Files", "*.xlsx")])
        if nombre_archivo:  # Si el usuario selecciona un archivo
            wb = openpyxl.load_workbook(nombre_archivo)
            sheet = wb.active
              # Suponiendo que esta es la lista donde almacenarás tus páginas

            i = 0
            for page in self.frames_notebook[:-1]:
                for j in range(len(page.float_vars)):
                    value = sheet.cell(row=i+1, column=2).value
                    # Si la valeur est présente dans la cellule
                    if value is not None:
                        # Met à jour la DoubleVar associée à l'Entry
                        page.float_vars[j].set(value)
                    i+=1

            print(f"Données chargées depuis : {nombre_archivo}")

    def set_text(self, entry, text):
        entry.delete(0,END)
        entry.insert(0,text)
   # def run_page(self):
    #   Run_program()
    


        
        

