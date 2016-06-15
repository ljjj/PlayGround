for i in range(332):
    for j in range((999-i)/2):
        a = i+1
        b = i+j+2
        c = 1000-a-b
        if a**2+b**2==c**2:
            print(str(a)+'*'+str(b)+'*'+str(c)+'='+str(a*b*c))
