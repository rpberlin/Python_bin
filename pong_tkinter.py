#!/usr/bin/env python3
import tkinter as tk
import random
import time

def moveBaseLR(base, dir, x, y = 0):
    x1, y1, x2, y2 = c.coords(base)
    if ((x1 > 0 and dir == 'l') or (x2 < 400 and dir == 'r')):
        c.move(base, x, y)
        c.update()

def startBall(ball, sp):
    s = random.randint(-sp, sp)
    x, y = s, 0-sp #Start. Ball to move in random direction. 0-sp is used to get negative value
    c.move(ball, x, y)

    for p in range(1, 500000):
        l, t, r, b = c.coords(ball)
        txtS.delete(0)
        txtS.insert(0, str(p))
        #Need to change direction on hitting wall. Eight options are there.
        if(r >= 400 and x >= 0 and y < 0): #Ball moving ↗ and hit right wall
            x, y = 0-sp, 0-sp
        elif(r >= 400 and x >= 0 and y >= 0): #Ball moving ↘ and hit right wall
            x, y = 0-sp, sp
        elif(l <= 0 and x < 0 and y < 0): #Ball moving ↖ and hit left wall
            x, y = sp, 0-sp
        elif(l <= 0 and x < 0 and y >= 0): #Ball moving ↙ and hit left wall
            x, y = sp, sp
        elif(t <= 0 and x > 0 and y < 0): #Ball moving ↗ and hit top wall
            x, y = sp, sp
        elif(t <= 0 and x < 0 and y < 0): #Ball moving ↖ and hit top wall
            x, y = 0-sp, sp
        elif(b >= 385): #Ball reached base level. Check if base touches ball
            tchPt = l + 10 #Size is 20. Half of it.
            bsl, bst, bsr, bsb = c.coords(base)
            if(tchPt >= bsl and tchPt <= bsr): #Ball touch base
                n = random.randint(-sp, sp)
                x, y = n, 0-sp
            else: #Hit bottom. Failed
                c.itemconfigure(lblID, state='normal')
                global isFailed
                isFailed = True
                break
        
        time.sleep(.025)
        c.move(ball, x, y)
        c.update()

def restart():
    global isFailed
    if(isFailed == True):
        isFailed = False
        c.itemconfigure(lblID, state='hidden')
        c.moveto(base, 150, 385)
        c.moveto(ball, 190, 365)
        startBall(ball, ballsp)



if __name__ == "__main__":
    root = tk.Tk()
    root.minsize(400,400)
    basesp = 10
    ballsp = 5
    global isFailed
    isFailed = False

    c = tk.Canvas(width=400, height=400, background='#a0aa00')
    c.pack()
    base = c.create_rectangle(150, 385, 250, 400, fill='blue', outline='blue')
    ball = c.create_oval(190, 365, 210, 385, fill='red', outline='red')    
    txtS = tk.Entry(c, text='0')
    txtScore = c.create_window(350, 0, anchor='nw', window=txtS)

    lblM = tk.Label(c, text='Failed!!!Press Enter key to start again')
    lblID = c.create_window(100, 190, anchor='nw', window=lblM)
    c.itemconfigure(lblID, state='hidden')

    #root.bind("<KeyPress-Left>", lambda event: moveBaseLR(base, 'l', 0-basesp))
    #root.bind("<KeyPress-Right>", lambda event: moveBaseLR(base, 'r', basesp))
    root.bind("<Left>", lambda event: moveBaseLR(base, 'l', 0-basesp))
    root.bind("<Right>", lambda event: moveBaseLR(base, 'r', basesp))
    root.bind("<Return>", lambda event: restart())
    
    startBall(ball, ballsp)

    root.mainloop()