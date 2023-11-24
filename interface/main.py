import tkinter as tk
from tkextrafont import Font
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap import Style
from isystem01 import isystem01
from isystem02 import isystem02
from isystem03 import isystem03
from isystem04 import isystem04
from isystem05 import isystem05
from main_page import main_page



class MainProgram(ttk.Frame):
    def __init__(self, master_window):
        super().__init__(master_window)
        self.pack(fill=BOTH, expand=YES)
        NavBar(self)
        DowBar(self)
        
        self.pages = [main_page(self),isystem01(self),isystem02(self),isystem03(self),isystem04(self),isystem05(self)]
        self.actual_page = 0
        self.show_pages(0)

    def show_pages(self, page):
        self.pages[self.actual_page].unshow()
        self.actual_page = page
        frame = self.pages[self.actual_page]
        frame.show()
    
    def show_values(self):
        t = self.pages[0].return_all_values()
        print(t)


class DowBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=BOTTOM, fill=BOTH)
        
        frame = ttk.Frame(self, style='secondary')
        frame.pack()
        
        sauvegarde = ttk.Button(frame, style='success.Solid.TButton', text='sauvegarder', width=35)
        sauvegarde.pack(side=LEFT)
        
        recuperer = ttk.Button(frame, style='danger.Solid.TButton', text='recuperer', width=35)
        recuperer.pack(side=LEFT)
        
        executer = ttk.Button(frame, style='primary.Solid.TButton', text='executer', width=35, command=lambda: frameprincipal.show_values())
        executer.pack(side=LEFT)



class NavBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=LEFT,fill=Y)
        
        system01 = ttk.Button(self, style='info.Solid.TButton', text='Main',width=15, command = lambda : frameprincipal.show_pages(0))
        system01.pack(side=TOP,padx=10,pady=[250,10])


        system01 = ttk.Button(self, style='info.Solid.TButton', text='PID fuelflow',width=15, command = lambda : frameprincipal.show_pages(1))
        system01.pack(side=TOP,padx=10,pady=[0,10])


        system02 = ttk.Button(self, style='info.Solid.TButton', text='Propeller',width=15, command = lambda : frameprincipal.show_pages(2))
        system02.pack(side=TOP,pady=[0,10])


        system03 = ttk.Button(self, style='info.Solid.TButton', text='P. Management',width=15, command = lambda : frameprincipal.show_pages(3))
        system03.pack(side=TOP,pady=[0,10])


        system04 = ttk.Button(self, style='info.Solid.TButton', text='Engine',width=15, command = lambda : frameprincipal.show_pages(4))
        system04.pack(side=TOP,pady=[0,10])

        
        system05 = ttk.Button(self, style='info.Solid.TButton', text='Engine MVEM',width=15, command = lambda : frameprincipal.show_pages(5))
        system05.pack(side=TOP)

        




if __name__ == "__main__":

    app = ttk.Window("PAPPL", "cyborg", resizable=(False, False), size=[1080,720], iconphoto='icon.png')
    app.place_window_center()

    font = Font(file="./interface/Nunito/static/Nunito-Black.ttf")

    style = Style()
    style.configure('main.TLabel', font=(font, 40))

    MainProgram(app)

    app.mainloop()
