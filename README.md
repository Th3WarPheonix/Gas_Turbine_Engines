# Gas_Turbine_Engines
Currently an amalgamation of gas turbine engine cycle design equations and component design equations from various sources

## Sources:  
[1] University class  
[2] Aircraft Propulsion 3e Farokhi  
[3] Aircraft Engines Design 2e Mattingly Heiser Pratt

## Legend
Listing of variables, subscripts, abbreviations, and symbols found in code and their basic descritpion for eqasy reference and consistency


### Variables 

| Variable Name | Actual Name | Description | Book Symbol
| --- | --- | --- | --- |
| wngld | wing loading | takeoff_weight/wing_planform_area | Wto/S
| thstld | thrust laoding | thrust/takeoff_weight | T/Wto
| engtype | engine type | tj, lbtf, tp, etc | 
| engmode | engine mode | wet or dry | 
| thrtlrto | throttle ratio |  | TR
| thstlps | installed thrust lapse | thrust/sea_level_thrust| thstlps
| instwf | instantaneous weight fraction | weight/takeoff_weight | instwf
|
| CLmax | maximum coefficient of lift | lift/density/planform_area/V^2/2 | C_Lmax
| CD | coefficient of drag | drag/density/planform_area/V^2/2 | C_D
| CD0 | parasite drag | drag coeff at zero lift <br> CDmin + K"Cmin^2| C_D0
| CDR | coefficient of additional drags | | C_DR
| CDstar | drag coefficient at max L/D | | C*_D
| Kp1 | invscid/induced drag coefficient | 1/pi/aspct_rto/efficiency | K'
| Kp2 | viscid drag coefficient | K" | K"
| K1 | drag polar coefficient 1 | K1 = K' + K" | K1
K2 | drag polar coefficient 2 | K2 = -2K"CLmin | K2
|
| drgcf | drag coefficient |  | &epsilon;
| frccf | friction coefficient |  | &mu;
| velrt | velocity ratio | velocity/stall_velocity | k
|
| thetat | nondimensional total temperature | Tt/Tstd | &theta;_0
| thetas | nondimensional static temperature | Ts/Tstd | &theta;
| deltat | nondimensional total pressure | Pt/Pstd | &delta;_0
| deltas | nondimensional static pressure | Ps/Pstd | &delta;
| sigmat | nondimensional total density | densityt/densitystd | &sigma;_0
| sigmat | nondimensional static density | densitys/densitystd | &psi;
| dynpress | dynamic pressure | density*velocity^2/2 | q
| thetabrk | theta break | control system maximum Tt4 and compressor pressure ratio | &theta;_0_break
| alt | altiude | | 
| gamma | ratio of specific heat | cp/cv | &gamma;
| Rgas | gas specific gas constant | | R
|
| swpang | wing sweep angle |  | &Lambda;
| ldfactor | load factor |  | n
| thstang | angle of thrust vector to chord line |  | &phi; 
| emptwf | empty aircraft weight fraction | We/Wto | &Gamma;
| drgthst | total drag-to-thrust ratio | | u
| tsfc | thrust specific fuel consumption | fuel_rate/thrust | tsfc
| mfp | mass flow parameter | fcn(gas, mach) | mfp 

### Abbreviations 
| Abbr | Expansion
| --- | --- |
| BCA | best cruise altitude
| BCM | best cruise altitude
| SLS | sea level static
| anlys | analysis

### Subscripts

| Sub | Name
| --- | --- |
| dry | afterburner off
| wet | afterburner on
| horz | horizontal
| accel | acceleration
| ln | landing
| to | takeoff
| td | touchdown
| obs | obstacle
| drgplr | drag polar
