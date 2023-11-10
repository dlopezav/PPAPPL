import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *


class isystem05(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')
        titre_text = "System 5"
        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.update_idletasks() 

    def unshow(self):
        self.pack_forget()