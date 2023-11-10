import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap import Style
from isystem01 import isystem01
from isystem02 import isystem02



class MainProgram(ttk.Frame):
    def __init__(self, master_window):
        super().__init__(master_window)
        self.pack(fill=BOTH, expand=YES)
        self.pages = [isystem01, isystem02]

        NavBar(self)
        
        titre_text = "Bienvenue à votre sistème"
        titre = ttk.Label(self, text=titre_text, style='main.TLabel')
        titre.pack(fill=X, pady=10, padx=20)

        titre_text = "Commencez par le système numero 1..."
        titre = ttk.Label(self, text=titre_text)
        titre.pack(fill=X, pady=10,padx=20)

    def show_pages(self,page):
        frame = self.pages[page]
        frame.tkraise()


class NavBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=LEFT,fill=Y)

        system01 = ttk.Button(self, style='info.Solid.TButton', text='Main',width=15)
        system01.pack(side=TOP,padx=10,pady=[300,10])


        system01 = ttk.Button(self, style='info.Solid.TButton', text='PID fuelflow',width=15, lambda : frameprincipal.show_pages() )
        system01.pack(side=TOP,padx=10,pady=[0,10])


        system02 = ttk.Button(self, style='info.Solid.TButton', text='Propeller',width=15)
        system02.pack(side=TOP,pady=[0,10])


        system03 = ttk.Button(self, style='info.Solid.TButton', text='P. Management',width=15)
        system03.pack(side=TOP,pady=[0,10])


        system04 = ttk.Button(self, style='info.Solid.TButton', text='Engine',width=15)
        system04.pack(side=TOP,pady=[0,10])

        
        system05 = ttk.Button(self, style='info.Solid.TButton', text='Engine MVEM',width=15)
        system05.pack(side=TOP)
        




if __name__ == "__main__":

    app = ttk.Window("PAPPL", "cyborg", resizable=(False, False), size=[1080,720], iconphoto='icon.png')
    app.place_window_center()

    style = Style()
    style.configure('main.TLabel', font=('Nunito Sans', 60))
    
    MainProgram(app)

    app.mainloop()
