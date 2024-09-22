USE DA_Dev
GO

create function nztm_to_lat( @N float, @E float) returns float
as

Begin 

Declare
-- Constants
@a      float = 6378137,
@f      float = 1/298.257222101,
@phizero float =    0,
@lambdazero float = 173,
@Nzero  float = 10000000,
@Ezero  float = 1600000,
@kzero  float = 0.9996,

-- Common variables
@b  float = 0,
@esq float  = 0,
@mzero float = 0,
@A0 float = 0,
@A2 float = 0,
@A4 float = 0,
@A6 float = 0,



-- Variables
@Nprime float = 0,
@mprime float = 0,
@smn    float = 0,
@G      float = 0,
@sigma  float = 0,
@phiprime   float = 0,
@rhoprime   float = 0,
@upsilonprime float = 0,
@psiprime   float = 0,
@tprime float = 0,
@Eprime float = 0,
@chi    float = 0,
@term_1 float = 0,
@term_2 float = 0,
@term_3 float = 0,
@term_4 float = 0,
@term1  float = 0,
@term2  float = 0,
@term3  float = 0,
@term4  float = 0,

--Output
@latitude   float = 0, 
@longitude  float = 0; 

--Calculation: To NZTM to NZGD2000/WGS84

--Common

SET @b          = @a*(1-@f)
SET @esq        = 2*@f-POWER(@f,2)
SET @A0         =1-@esq/4-3*POWER(@esq,2)/64-5*POWER(@esq,3)/256 
SET @A2         =0.375*(@esq+POWER(@esq,2)/4+15*power(@esq,3)/128) 
SET @A4         =15*(POWER(@esq,2)+3*POWER(@esq,3)/4)/256  
SET @A6         =35*POWER(@esq,3)/3072 

SET @Nprime     = @N-@Nzero
SET @mprime     = @Nprime/@kzero
SET @smn        = (@a-@b)/(@a+@b)
SET @G          = @a*(1-@smn)*(1-POWER(@smn,2))*(1+9*POWER(@smn,2)/4+225*POWER(@smn,4)/64)*PI()/180.0
SET @sigma      = @mprime*PI()/(180*@G)
SET @phiprime   = @sigma+(3*@smn/2-27*POWER(@smn,3)/32)*SIN(2*@sigma)+(21*POWER(@smn,2)/16-55*POWER(@smn,4)/32)*SIN(4*@sigma)+(151*POWER(@smn,3)/96)*SIN(6*@sigma)+(1097*POWER(@smn,4)/512)*SIN(8*@sigma)
SET @rhoprime   = @a*(1-@esq)/POWER((1-@esq*POWER((SIN(@phiprime)),2)),1.5)
SET @upsilonprime   =@a/SQRT(1-@esq*POWER((SIN(@phiprime)),2)) 

SET @psiprime   = @upsilonprime/@rhoprime
SET @tprime     = TAN(@phiprime)
SET @Eprime     = @E-@Ezero
SET @chi        = @Eprime/(@kzero*@upsilonprime)
SET @term_1     = @tprime*@Eprime*@chi/(@kzero*@rhoprime*2)
SET @term_2     = @term_1*POWER(@chi,2)/12*(-4*POWER(@psiprime,2)+9*@psiprime*(1-POWER(@tprime,2))+12*POWER(@tprime,2))
SET @term_3     = @tprime*@Eprime*POWER(@chi,5)/(@kzero*@rhoprime*720)*(8*POWER(@psiprime,4)*(11-24*POWER(@tprime,2))-12*POWER(@psiprime,3)*(21-71*POWER(@tprime,2))+15*POWER(@psiprime,2)*(15-98*POWER(@tprime,2)+15*POWER(@tprime,4))+180*@psiprime*(5*POWER(@tprime,2)-3*POWER(@tprime,4))+360*POWER(@tprime,4))
SET @term_4     = @tprime*@Eprime*POWER(@chi,7)/(@kzero*@rhoprime*40320)*(1385+3633*POWER(@tprime,2)+4095*POWER(@tprime,4)+1575*POWER(@tprime,6))
SET @term1      = @chi*(1/COS(@phiprime))
SET @term2      = POWER(@chi,3)*(1/COS(@phiprime))/6*(@psiprime+2*POWER(@tprime,2))
SET @term3      = POWER(@chi,5)*(1/COS(@phiprime))/120*(-4*POWER(@psiprime,3)*(1-6*POWER(@tprime,2))+POWER(@psiprime,2)*(9-68*POWER(@tprime,2))+72*@psiprime*POWER(@tprime,2)+24*POWER(@tprime,4))
SET @term4      = POWER(@chi,7)*(1/COS(@phiprime))/5040*(61+662*POWER(@tprime,22)+1320*POWER(@tprime,4)+720*POWER(@tprime,6))

SET @latitude   = (@phiprime-@term_1+@term_2-@term_3+@term_4)*180/PI()
SET @longitude  = @lambdazero+180/PI()*(@term1-@term2+@term3-@term4)



    return @latitude;
end
go

create function nztm_to_lONG( @N float, @E float) returns float
as

Begin 

Declare
-- Constants
@a      float = 6378137,
@f      float = 1/298.257222101,
@phizero float =    0,
@lambdazero float = 173,
@Nzero  float = 10000000,
@Ezero  float = 1600000,
@kzero  float = 0.9996,

-- Common variables
@b  float = 0,
@esq float  = 0,
@mzero float = 0,
@A0 float = 0,
@A2 float = 0,
@A4 float = 0,
@A6 float = 0,



-- Variables
@Nprime float = 0,
@mprime float = 0,
@smn    float = 0,
@G      float = 0,
@sigma  float = 0,
@phiprime   float = 0,
@rhoprime   float = 0,
@upsilonprime float = 0,
@psiprime   float = 0,
@tprime float = 0,
@Eprime float = 0,
@chi    float = 0,
@term_1 float = 0,
@term_2 float = 0,
@term_3 float = 0,
@term_4 float = 0,
@term1  float = 0,
@term2  float = 0,
@term3  float = 0,
@term4  float = 0,

--Output
@latitude   float = 0, 
@longitude  float = 0; 

--Calculation: To NZTM to NZGD2000/WGS84

--Common

SET @b          = @a*(1-@f)
SET @esq        = 2*@f-POWER(@f,2)
SET @A0         =1-@esq/4-3*POWER(@esq,2)/64-5*POWER(@esq,3)/256 
SET @A2         =0.375*(@esq+POWER(@esq,2)/4+15*power(@esq,3)/128) 
SET @A4         =15*(POWER(@esq,2)+3*POWER(@esq,3)/4)/256  
SET @A6         =35*POWER(@esq,3)/3072 

SET @Nprime     = @N-@Nzero
SET @mprime     = @Nprime/@kzero
SET @smn        = (@a-@b)/(@a+@b)
SET @G          = @a*(1-@smn)*(1-POWER(@smn,2))*(1+9*POWER(@smn,2)/4+225*POWER(@smn,4)/64)*PI()/180.0
SET @sigma      = @mprime*PI()/(180*@G)
SET @phiprime   = @sigma+(3*@smn/2-27*POWER(@smn,3)/32)*SIN(2*@sigma)+(21*POWER(@smn,2)/16-55*POWER(@smn,4)/32)*SIN(4*@sigma)+(151*POWER(@smn,3)/96)*SIN(6*@sigma)+(1097*POWER(@smn,4)/512)*SIN(8*@sigma)
SET @rhoprime   = @a*(1-@esq)/POWER((1-@esq*POWER((SIN(@phiprime)),2)),1.5)
SET @upsilonprime   =@a/SQRT(1-@esq*POWER((SIN(@phiprime)),2)) 

SET @psiprime   = @upsilonprime/@rhoprime
SET @tprime     = TAN(@phiprime)
SET @Eprime     = @E-@Ezero
SET @chi        = @Eprime/(@kzero*@upsilonprime)
SET @term_1     = @tprime*@Eprime*@chi/(@kzero*@rhoprime*2)
SET @term_2     = @term_1*POWER(@chi,2)/12*(-4*POWER(@psiprime,2)+9*@psiprime*(1-POWER(@tprime,2))+12*POWER(@tprime,2))
SET @term_3     = @tprime*@Eprime*POWER(@chi,5)/(@kzero*@rhoprime*720)*(8*POWER(@psiprime,4)*(11-24*POWER(@tprime,2))-12*POWER(@psiprime,3)*(21-71*POWER(@tprime,2))+15*POWER(@psiprime,2)*(15-98*POWER(@tprime,2)+15*POWER(@tprime,4))+180*@psiprime*(5*POWER(@tprime,2)-3*POWER(@tprime,4))+360*POWER(@tprime,4))
SET @term_4     = @tprime*@Eprime*POWER(@chi,7)/(@kzero*@rhoprime*40320)*(1385+3633*POWER(@tprime,2)+4095*POWER(@tprime,4)+1575*POWER(@tprime,6))
SET @term1      = @chi*(1/COS(@phiprime))
SET @term2      = POWER(@chi,3)*(1/COS(@phiprime))/6*(@psiprime+2*POWER(@tprime,2))
SET @term3      = POWER(@chi,5)*(1/COS(@phiprime))/120*(-4*POWER(@psiprime,3)*(1-6*POWER(@tprime,2))+POWER(@psiprime,2)*(9-68*POWER(@tprime,2))+72*@psiprime*POWER(@tprime,2)+24*POWER(@tprime,4))
SET @term4      = POWER(@chi,7)*(1/COS(@phiprime))/5040*(61+662*POWER(@tprime,22)+1320*POWER(@tprime,4)+720*POWER(@tprime,6))

SET @latitude   = (@phiprime-@term_1+@term_2-@term_3+@term_4)*180/PI()
SET @longitude  = @lambdazero+180/PI()*(@term1-@term2+@term3-@term4)



    return @lONGITUDE;
end
go