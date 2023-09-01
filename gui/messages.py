import os

class TerminalMessages:
    """ Messages that print in terminal """
    def __init__(self):
        tsize = os.get_terminal_size()
        self.xsize = tsize.columns

    def div_line(self):
        pline = " " + '*'*(self.xsize - 2) + " "
        print(pline)

    def centered_text(self, message=None):
        if message is not None:
            n_space = self.xsize - 8 - len(message)
            if (n_space % 2) == 0:
                n_left = int(n_space/2)
                n_right = int(n_space/2)
            else:
                n_left = int(n_space/2)
                n_right = int(n_space/2+1)
            pline =  " ***" + " "*n_left + message +" "*n_right + "*** "
        else:
            n_space = self.xsize - 8
            pline =  " ***" + " "*n_space + "*** "
        print(pline)

    def boxed_message(self, message):
        print('\n')
        self.div_line()
        self.div_line()
        self.centered_text()
        self.centered_text(message)
        self.centered_text()
        self.div_line()
        self.div_line()
        print('\n')

    def welcome(self):
        self.boxed_message("Welcome to the Dysh GUI")

    def goodbye(self):
        self.boxed_message("Goodbye!")