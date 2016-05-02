ptf @
* parameter data
s21   = @s21      @   
s22   = @s22      @
s24    =@s24      @
s26    =@s26      @
a21    =@a21      @
a22    =@a22      @
a24    =@a24      @
a26    =@a26      @
c11    =@c11     @
c12    =@c12     @
c14    =@c14     @
c16    =@c16     @
lhg    =@lhg       @
chc    =@chc @
nmac    =@   nmac  @  
mac_r  =@  mac_r @
LKs    =@LKs      @
Koc    =@Koc  @ 
half   =@half@
h2    =@h2 @
expr   =@expr @
EST    =@EST @
rsw    =@rsw   @
rst    =@rst  @
Ea     =@Ea     @
wsm    =@Wsm  @
mrad   =@mrad   @
h_3    =@h_3 @

h_2=half*h2
mac=0.5*nmac*3.14159*mac_r^2

c13=c12*1.40496
c15=c14*1.26087
c17=c16*1.0

s23=s22
s25=s24
s27=s26

a23=a22
a25=a24
a27=a26

s11  = s21
s12  = s22  
s13  = s23  
s14  = s24  
s15=s25
s16=s26
s17=s27

n21  = 2+3*a21
n22  = 2+3*a22
n23  = 2+3*a23
n24  = 2+3*a24
n25=2+3*a25
n26=2+3*a26
n27=2+3*a27

c21  = c11*s11^n21
c22  = c12*s12^n22
c23  = c13*s13^n23
c24  = c14*s14^n24
c25=c15*s15^n25
c26=c16*s16^n26
c27=c17*s17^n27

* template and model input data
C:\PEST\nash_pest\rzdat_deg.tpl    C:\RZWQM2\Nasha_Fox\nashua13\rzwqm.dat
