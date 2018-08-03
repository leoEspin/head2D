#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Testing GUIs and OOP with python
Created on Thu Aug  2 14:41:24 2018

@author: Leo
"""
import sys
from PyQt5.QtWidgets import (QWidget, QApplication, QLabel, QPushButton, 
                             QLineEdit)
from PyQt5.QtCore import QSize 

class InputFile(QWidget):
    
    def __init__(self):
        super().__init__()
        
        self.initUI()
          
    def initUI(self):      
        self.lbl1 = QLabel(self)  
        self.lbl1.setText('Points per unit length (X)')
        lbl2 = QLabel('Points per unit length (Y)', self)        
        lbl3 = QLabel('Points per time unit', self)       
        lbl4 = QLabel('Diffusion coefficient', self)        
        lbl5 = QLabel('Final time', self)    
        self.le1 = QLineEdit('100',self)
        self.le2 = QLineEdit('100',self)
        self.le3 = QLineEdit('1e5',self)
        self.le4 = QLineEdit('1.0',self)
        self.le5 = QLineEdit('1.0',self)

        self.lbl1.move(15, 10)
        self.le1.move(180, 10)     
        lbl2.move(15, 40)
        self.le2.move(180, 40)
        lbl3.move(15, 70)
        self.le3.move(180, 70)
        lbl4.move(15, 100)
        self.le4.move(180, 100)
        lbl5.move(15, 130)
        self.le5.move(180, 130)

        self.btn1 = QPushButton('Ok', self)
        self.btn1.move(35, 170)
        self.btn1.clicked.connect(self.clickMethod)
        
        self.btn2 = QPushButton('Cancel', self)
        self.btn2.move(200, 170)     
        self.btn2.clicked.connect(self.close)
        
        self.setGeometry(300, 300, 320, 245)
#        self.setMinimumSize(QSize(320, 245))   
        self.setFixedSize(320, 245)
        self.setWindowTitle('Input parameters')
        
    def clickMethod(self):
        file = open('input','w')
        file.write(self.le1.text() + '\n')
        file.write(self.le2.text() + '\n')
        tmp=float(self.le3.text())
        tmp=1./tmp
        tmp=str(tmp)
        file.write(tmp + '\n')
        file.write(self.le4.text() + '\n')
        file.write(self.le5.text() + '\n')
        file.close
        sys.exit()
               
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    objetito = InputFile()
    objetito.show()
    sys.exit(app.exec_())

