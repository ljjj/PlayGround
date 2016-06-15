def palindrome(s):
    p = True
    for k in range(len(s)/2):
        p = p and s[k]==s[len(s)-1-k]
        if not p:
            break
    return p

largest = 0
for i in range(900):
    for j in range(900-i):
        x = 999-i
        y = 999-i-j
        z= (999-i)*(999-i-j)
        if palindrome(str(z)):
            print(str(x)+'*'+str(y)+'='+str(z))
            if z > largest:
                largest = z
