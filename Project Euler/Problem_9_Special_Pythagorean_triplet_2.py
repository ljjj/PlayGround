N = 1000
for i in range(N/3):
    a = i+1
    b = (N**2.0-2*N*a)/(2*N-2*a)
    if b.is_integer():
        b = int(b)
        c = N-a-b
        print(str(a)+'*'+str(b)+'*'+str(c)+'='+str(a*b*c))
