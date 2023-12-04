sqrtLim = 500
sqrt2 = 1.4142135623730950488 

def sqrt(n):
    #approx
    odd = 0
    a = int(n).bit_length()-1
    if a & 1: odd = 1
    a >>= 1 #div by 2 / cuts out odd
    lw = 2**a #lowerbound
    up = lw* [sqrt2, 2][odd] #syntax shenanigans

    for _ in range(sqrtLim):
        m = (up+lw)/2
        t = m**2 -n
        if t == 0: return(m)
        elif t < 0: lw = m
        else: up = m
    return(m)

loopLim = 1000

def sgn(n):
    if n > 0: return(1)
    elif n < 0: return(-1)
    return(0)

class polynom:
    def __init__(self, cef, NaN = False):
        self.cef = cef
        self.NaN = NaN
    
    @property
    def order(self):
        return( len(self.cef)-1 )
    
    def __call__(self, val):
        s = 0
        for i in range(self.order+1):
            s += self[i] * (val**i)
        return(s)

    def __str__(self):
        if self.NaN: return("NaN")

        txt = ""
        for i, c in enumerate(self.cef):
            #zero case val
            if c == 0: continue
            #sign
            elif c > 0 and i != 0: txt += " + "
            elif c < 0 and i != 0: txt += " - " ; c = str(c)[1:]
            #zero/one case exp
            if i not in (self.order, self.order-1):
                ord = self.order - i
                if ord < 10: txt += f"{c}x^{ord}"
                else: txt += f"{c}x^"+"{"+str(ord)+"}"
            elif i == self.order-1: txt += f"{c}x"
            elif i == self.order: txt += f"{c}"

        return(txt)
    
    def __getitem__(self, indx):
        if indx > self.order: raise(Exception("Faulty index"))
        return(self.cef[self.order - indx])
    
    def __iter__(self):
        return(iter(self.cef))

def ddx(p : polynom): 
        if p.NaN: raise(Exception("Cannot differiantate NaNs"))
        #const case
        if p.order == 0: return(0)

        Dcefs = []
        for i in range(0, p.order):
            i = p.order-i
            Dcefs.append(i*p[i])
        
        Rddx = polynom(Dcefs)
        return(Rddx)

def quadratic(p2 : polynom):
    if p2.order != 2: raise(Exception("Quadratic only handles polynomials of 2nd degree"))
    d = p2[1]**2 - 4*p2[2]*p2[0]
    if d < 0: return( [] )
    
    elif d == 0:
        return( [-p2[1]/(2*p2[2])] )

    d = sqrt(d)
    ans = []
    for s in (1, -1):
        c = -p2[1]+s*d
        c /= 2*p2[2]
        ans.append(c)

    ans = sorted(ans)

    return(ans)

def newtonIter(p, g):

    fp = ddx(p)
    for _ in range(loopLim):
        g -= p(g)/fp(g)

    return(g)

#fix this shit
def bisectSolve(p, ls): 
    a, b = sorted(ls)
    if sgn( p(a) ) == sgn( p(b) ): raise(Exception("Cannot bisect from section sharing sign."))
    if p(b) < p(a):
        a, b = b, a
    
    for _ in range(loopLim):

            m = (a+b)/2
            tmp = p(m)
            if tmp == 0: return(m)
            elif tmp < 0: a = m
            elif tmp > 0: b = m
    
    return(m)

def findCP(p):
    df = ddx(p)
    CP = findRoot(df)
    CP = sorted(CP)
    return(CP)
    

def findRoot(p : polynom):
    if p.order == 2: return( quadratic(p) )
    elif p.order == 1: return(-p[0]/p[1])
    
    Zs = [] #here go the zeros

    CP = findCP(p)
    Vs = []; Qs = []
    ddf = ddx(ddx(p))

    for P in CP:
        Q = sgn( ddf(P) ) 
        if Q == 0: CP.remove(P) ; break #terassi piste check

        Qs.append(Q)
        Vs.append( p(P) )
    

    if len(CP) == 0:
        #ADD CHECKS AND FAIL SAFES, HERE AND IN NEWTONITER
        m = newtonIter(p, 1)
        Zs.append(m)
        return(Zs)
    
    #edge checks
    for indx in (0, -1):
        cs = sgn(Vs[indx])
        if cs == Qs[indx]: continue
        
        if indx == 0: s = -1
        else: s = 1

        mv = 1
        bnd = CP[indx]+s*mv
        while cs == sgn(p(bnd)): 
            mv <<= 1
            bnd = CP[indx]+s*mv

        m = bisectSolve(p, (bnd, CP[indx]))
        Zs.append(m)
    
    for indx in range(len(CP)-1):
        if sgn(Vs[indx]) == sgn(Vs[indx+1]): continue
        m = bisectSolve(p, (CP[indx], CP[indx+1]))
        Zs.append(m)

    Zs = sorted(Zs)
    return(Zs)


t = polynom((3, 2,1,-100, 3, 4, 0, 1, -3, 6, 0, 0, -5, -1, 3))
print(t)
print("Root(s):", findRoot(t))
