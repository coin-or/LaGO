#  MINLP written by GAMS Convert at 10/02/06 16:20:40
#  
#  Equation counts
#      Total        E        G        L        N        X        C
#          6        3        0        3        0        0        0
#  
#  Variable counts
#                   x        b        i      s1s      s2s       sc       si
#      Total     cont   binary  integer     sos1     sos2    scont     sint
#          6        3        3        0        0        0        0        0
#  FX      0        0        0        0        0        0        0        0
#  
#  Nonzero counts
#      Total    const       NL      DLL
#         17       15        2        0
# 
#  Reformulation has removed 1 variable and 1 equation


var x1 := 1, >= 0, <= 10;
var x2 := 1, >= 0, <= 10;
var b3 binary >= 0, <= 1;
var b4 binary >= 0, <= 1;
var b5 binary >= 0, <= 1;

minimize obj:    2*x1 + 3*x2 + 1.5*b3 + 2*b4 - 0.5*b5;

subject to

e2: x1^2 + b3 = 1.25;

e3: x2**1.5 + 1.5*b4 = 3;

e4:    x1 + b3 <= 1.6;

e5:    1.333*x2 + b4 <= 3;

e6:  - b3 - b4 + b5 <= 0;
