# Import necessary libraries
import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd
import math
from Notebook import *

# Define the main_page class
class main_wthout(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')

        # Title for the main page
        titre_text = "MAIN WITH WTHOUT"
        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

        # Initialize values, texts, and notebook frames
        self.values = []
        self.texts = []
        
        # Add notebook frames to the notebook
        for count, text in enumerate(self.texts):
            self.notebook.add(self.frames_notebook[count], text=text)

        # Create a frame with buttons for saving, retrieving, and executing
        self.frame = ttk.Frame(self, style='secondary')
        self.sauvegarde = ttk.Button(self.frame, command=lambda: self.sauvegarder_dans_excel(), style='success.Solid.TButton', text='sauvegarder', width=35) 
        self.recuperer = ttk.Button(self.frame, command=lambda: self.charger_depuis_excel(), style='danger.Solid.TButton', text='recuperer', width=35)
        self.executer = ttk.Button(self.frame, style='primary.Solid.TButton', text='executer', width=35, command=lambda: frameprincipal.executer())
        
        # Pack the buttons
        self.frame.pack(side=BOTTOM)
        self.recuperer.pack(side=LEFT)
        self.sauvegarde.pack(side=LEFT)
        self.executer.pack(side=LEFT)

    # Show the main page
    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.update_idletasks()

    # Hide the main page
    def unshow(self):
        self.pack_forget()

    # Save data to an Excel file - second column with all the values
    def sauvegarder_dans_excel(self):
        return
       

    # Load data from an Excel file - data must be saved by "sauvegarder_dans_excel"
    def charger_depuis_excel(self):
       return

    # Set text for an entry
    def set_text(self, entry, text):
        entry.delete(0, END)
        entry.insert(0, text)

    def executer_program(self):
        
        import interface.interface as run_program
        run_program


    


        