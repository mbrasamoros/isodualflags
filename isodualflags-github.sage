# This is a set of sage commands to deal with isometry-dual flags of codes
# Sage is a free open-source mathematics software system licensed under the GPL. It can be downloaded from www.sagemath.org, or used online through the link https://sagecell.sagemath.org/
# The commands here presented have been used to compute examples for the references
# [1] Bras-Amorós, Maria; Castellanos, Alonso Sepúlveda; Quoos, Luciane: The isometry-dual property in flags of two-point algebraic geometry codes. IEEE Trans. Inform. Theory 68 (2022), no. 2, 828–838. 
# [2] Bras-Amorós, Maria; Castellanos, Alonso Sepúlveda; Quoos, Luciane: Isometry-dual Flags of many-point AG codes, submitted (2022). 


print("##############################")
print("LOADING FUNCTIONS")

def Hb(n,g,D,P,Q,QQ,b,bb,mina,maxa):
  olddim=0
  Hbset=set({})
  ini=-b-bb
  if mina>ini:
    ini=mina
  fin=2*g-b-bb
  if fin > maxa:
    fin=maxa
  for a in [ini..fin]:
    newdim= (a*P+b*Q+bb*QQ).dimension()
    if newdim != olddim:
      Hbset = Hbset.union({a})
      olddim=newdim
  for a in [fin+1..maxa]:
    Hbset = Hbset.union({a})
  return Hbset

def Hstar(n,g,D,P,Q,QQ,b,bb,mina,maxa):
  HH=sorted(Hb(n,g,D,P,Q,QQ,b,bb,mina,n+2*g-1-b-bb))
  if len(HH) >= 1:
    Hbstarset=set({HH[0]})
    HH=[HH[i] for i in [1..len(HH)-1]]
    olddim=1
  else:
    Hbstarset=set({})
  for a in HH:
    C=codes.EvaluationAGCode(D, a*P+b*Q+bb*QQ)
    newdim=C.dimension()
    if newdim != olddim:
      Hbstarset = Hbstarset.union({a})
      olddim=newdim
  return Hbstarset

def fastfindx(FF,n,g,D,P,Q,QQ,b,bb):
  Hpreaux=sorted(Hstar(n,g,D,P,Q,QQ,b,bb,-b-bb,n+2*g-1-b-bb))
  Haux=[Hpreaux[i] for i in [0..len(Hpreaux)-2] ];
  l=len(Haux)
  if l == 0:
    print ("empty code flag for b=",b)
    return False, V(0)
  V=VectorSpace(FF,n)
  rows=[]
  for i in [0..floor(l/2)-1]:
    primmatrix=codes.EvaluationAGCode(D,Haux[i]*P+b*Q+bb*QQ).generator_matrix()
    dualmatrix=codes.EvaluationAGCode(D,Haux[l-i-1]*P+b*Q+bb*QQ).generator_matrix()
    temprows=flatten([V([primmatrix[ii][kk]*dualmatrix[jj][kk] for kk in [0..n-1]]) for ii in [0..i] for jj in [0..l-i-1]])
    rows=rows+temprows
    N=(matrix(rows)).transpose().kernel()
    if N.dimension() == 0:
      return False, V(0)
    elif N.dimension() == 1:
      x=N.basis()[0]
      if FF(0) not in x:
        return True, x
      else:
        return False, V(0)
  V=VectorSpace(FF,N.dimension())
  print(V)
  for v in V:
    x=v*N.basis_matrix()
    if (FF(0) not in x):
      return True, x
    return False, V(0)


print("##############################")
print("HERMITE EXAMPLE (y^3=x^2+x)")

FF.<alpha>=GF(4)
K.<x> = FunctionField(FF); R.<t> = K[]
L.<y> = K.extension(t^3 -x^2-x)
P= L.places_infinite()[0]
places=L.places_finite()
n=len(places)-2
Q= places[0]
QQ=places[n+1]
D=[places[i] for i in [1..n]]
g=1
b=-1
bb=2
print("Fix the parameters")
print("D=",[[x.evaluate(D[i]),y.evaluate(D[i])] for i in [0..(n-1)]])
print("P=",P)
print("Q=",Q)
print("QQ=",QQ)
print("b=",b)
print("bb=",bb)


print("The Hb set corresponding to the flag of codes C(D,aP+bQ+bbQQ) is ")
print(Hb(n,g,D,P,Q,QQ,b,bb,-b-bb,n+2*g-1-b-bb))
print("The Hb^ set corresponding to the flag of codes C(D,aP+bQ+bbQQ) is ")
print(Hstar(n,g,D,P,Q,QQ,b,bb,b-bb,n+2*g-1-b-bb))

print("Check whether the flag of codes C(D,aP+bQ+bbQQ) satisfies the isometry-dual property and in case it satisfies it, give the corresponding isometry vector")
print(fastfindx(FF,n,g,D,P,Q,QQ,b,bb))
