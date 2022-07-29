#include atoms

#DEFVAR
    DUMM1       = IGNORE;
    DUMM2       = IGNORE;
    CO2		= C + 2O ;
    H2		= 2H ;                          {molecular hydrogen}
    NO		= N + O ;			{nitric oxide}
    NO2		= N + 2O ;			{nitrogen dioxide}
    HNO3	= H + N + 3O ;			{nitric acid}
    NH3		= N + 3H ;			{ammonia}
    SO2		= S + 2O ;			{sulfur dioxide}
    SO3		= S + 3O ;			{sulfur trioxide}
    HOSO2       = S + 3O + H ;
    H2SO4	= S + 4O + 2H ;
    O3		= 3O ;				{ozone}
    CH4		= C + 4H ;			{methane}
    C2H6	= 2C + 6H ;			{ethane}
    C3H8	= 3C + 8H ;			{propane}
    ALKA	= 			IGNORE; {>C2 alkanes}
    ETHE	= 			IGNORE; {ethene}
    ALKE	= 			IGNORE; {>C2 alkenes}
    AROM	= 			IGNORE;
    ACO2	= 2H + C + 2O ;			{ HCOOH }
    ACTA	= 2C + 4H + 2O ; 		{ CH3COOH }
    HCHO	= 2H + C + O ;			{formaldehyde}
    ALD2	= 			IGNORE;
    CH3OH       = C + 4H + O ;                  {methanol}
    C2H5OH      = 2C + 6H + O ;                 {ethanol}
    H2O2	= 2H + 2O ;			{hydrogen peroxide}
    ROOH	= 			IGNORE;
    HONO	= H + 2O + N ;			{nitrous acid}
    PAN		= 2C + 3H + 5O + N ;		{ CH3CO3NO2 }
    TPAN	= 4C + 3H + 6O + N ;		{ CHOCH=CHCO3NO2 }
    KET		= 			IGNORE;
    CRES	= 			IGNORE;
    DIAL	= 			IGNORE;
    GLYX	= 2C + 2H + 2O ;		{ CHOCHO }
    MGLY	= 3C + 4H + 2O ; 	 	{ CH3COCHO }
    NH4NO3	= 2N + 4H + 3O ;         	{ or AMNT }
    HCl		= H + Cl ;			{hydrogen chloride}
    R3N2	= N + IGNORE ;
    RAN1	= N + IGNORE ;
    RAN2	= N + IGNORE ;
    N2O5	= 2N + 5O ;			{dinitrogen pentoxide}
    HNO4	= H + N + 4O ;			{pernitric acid}
    NO3		= N + 3O ;	                {nitrogen trioxide}
    DMS		= 2C + 6H + S ; 		{ CH3-S-CH3 }
    DMSO	= 2C + 6H + S + O ;		{DMSO: CH3SOCH3}
    DMSO2       = 2C + 6H + S + 2O ;	        {DMSO2}
    DMOO        = 2C + 5H + S + 2O;             {CH3SCH2OO}   
    CH3S        = C + 3H + S;                   {}
    CH3SO       = C + 3H + S + O;               {}
    CH3SO2	= C + 3H + S + 2O ;		{NOT CH3SOO}
    CH3SO3	= C + 3H + S + 3O ;		{}
    CH3SO2H     = C + 4H + S + 2O ;		{CH3S(O)OH, MSIA}
    CH3SO3H	= C + 4H + S + 3O ;		{CH3S(OO)OH, MSA}
    CO		= C + O;			{carbon monoxide}
    HOCl	= H + O + Cl ;			{hydrochlorous acid}
    ClONO       = Cl + 2O + N ;                 {}
    ClNO2	= Cl + 2O + N ;			{chlorine nitrite}
    ClNO3	= Cl + 3O + N ;			{chlorine nitrate}
    Cl2		= 2Cl ;				{molecular chlorine}
    ClO3        = Cl + 3O ;                     {}
    Cl2O2       = 2Cl + 2O ;                    {}
    Cl2O3       = 2Cl + 3O ;                    {}
    HBr		= H + Br ;			{hydrogen bromide}
    HOBr	= H + O + Br ;			{hydrobromous acid}
    BrNO2	= Br + N + 2O ;			{bromine nitrite}
    BrNO3	= Br + N + 3O ;			{bromine nitrate}
    Br2		= 2Br ;				{molecular bromine}
    Br2O	= 2Br + O ;			{}
    BrCl	= Br + Cl ;			{bromine chloride}
    HI		= H + I ;			{hydrogen iodide}
    HOI		= H + O + I ;			{hypoiodous acid}
    HIO2	= H + I + 2O ;			{jjb was missing; only produced in the gas phase, warning! To be investigated}
    I2O2	= 2O + 2I ;			{}
    INO         = I + N + O ;                   {}
    INO2	= I + N + 2O ;			{iodine nitrite}
    INO3	= I + N + 3O ;			{iodine nitrate}
    I2		= 2I ;				{molecular iodine}
    ICl		= I + Cl ;			{iodine chloride}
    IBr		= I + Br ;			{iodine bromide}
    CH3I	= C + 3H + I ;			{iodomethane}
    CH2I2	= C + 2H + 2I ;			{diiodomethane}
    CH2ClI	= C + 2H + Cl + I ;		{chloroiodomethane}
    C3H7I	= 3C + 7H + I ;			{2-iodopropane} 
    CH2BrI      = C + 2H + Br + I ;             {bromoiodomethane}
    CHBr2I      = C + H + 2Br + I ;             {dibromoiodomethane}
    C2H5I       = 2C + 5H + I ;                 {iodoethane}
    HIO3        = H + I + 3O ;                  {}
    I2O         = I + 2O ;                      {}
    I2O3        = 2I + 3O ;                     {}
    I2O4        = 2I + 4O ;                     {}
    I2O5        = 2I + 5O ;                     {}

    NHS         =                       IGNORE; {generic species containing N and/or S}
    SOR         =                       IGNORE; {generic sulphur-containing oxygenated}
    SPAN        =                       IGNORE; {generic sulphur-containing PAN}
    RCl         =                       IGNORE; {generic chlorine-containing hydrocarbon}
    RBr         =                       IGNORE; {generic bromine-containing hydrocarbon}
    XOR         =                       IGNORE; {generic halogen-containing oxygenated}

{    PRN2	= 2N + IGNORE ;}
{    PRPN	=  N + IGNORE ;}
{    OZID	=      IGNORE ;}


{liquid phase}
 {_l1: small aerosol}
    NOl1	= IGNORE;			{nitric oxide}
    NO2l1	= IGNORE;			{nitrogen dioxide}
    HNO3l1	= IGNORE;			{nitric acid}
    NH3l1	= IGNORE;			{N + 3Hammonia}
    SO2l1	= IGNORE;			{sulfur dioxide}
    SO4l1	= IGNORE;  
    O3l1	= IGNORE;			{ozone}

    ACO2l1	= IGNORE; 			{ HCOOH }
    ACTAl1	= IGNORE;             		{ CH3COOH }
    HCHOl1	= IGNORE;			{formaldehyde}
    ALD2l1	= 			IGNORE;
    H2O2l1	= IGNORE;		{hydrogen peroxide}
    ROOHl1	= 			IGNORE;   
    HONOl1	= IGNORE;		{nitrous acid}

    HCll1	= IGNORE;			{hydrogen chloride}

    N2O5l1	= IGNORE;		{dinitrogen pentoxide}    
    HNO4l1	= IGNORE;			{pernitric acid}
    NO3l1	= IGNORE;		{nitrogen trioxide}                            
    DMSl1	= IGNORE;			{ CH3-S-CH3 }   
    DMSOl1	= IGNORE; 		{DMSO: CH3SOCH3}
    CH3SO2l1	= IGNORE; 		{}
    CH3SO3l1	= IGNORE; 		{}
    CH3SO3Hl1	= IGNORE; 		{}
    DMOOl1	= IGNORE;		{CH3SCH2OO}
    CH3Sl1	= IGNORE;		{}
    CH3SOl1	= IGNORE;		{}
    CH3SO2Hl1	= IGNORE;		{MSIA}
    DMSO2l1	= IGNORE;		{}

    COl1	= IGNORE; 			{carbon monoxide}
    HOCll1	= IGNORE;			{hydrochlorous acid}
    ClNO2l1	= IGNORE;			{chlorine nitrite}
    ClNO3l1	= IGNORE;			{chlorine nitrate}
    Cl2l1	= IGNORE;				{molecular chlorine}
    HBrl1	= IGNORE; 			{hydrogen bromide}
    HOBrl1	= IGNORE;			{hydrobromous acid}
    BrNO2l1	= IGNORE; 			{bromine nitrite}
    BrNO3l1	= IGNORE; 			{bromine nitrate}
    Br2l1	= IGNORE; 				{molecular bromine}
    BrCll1	= IGNORE; 			{bromine chloride}
    HIl1	= IGNORE; 			{hydrogen iodide}
    HOIl1	= IGNORE; 			{hypoiodous acid}
    I2O2l1	= IGNORE; 			{}
    INO2l1	= IGNORE; 			{iodine nitrite}
    INO3l1	= IGNORE; 		{iodine nitrate}
    I2l1	= IGNORE; 				{molecular iodine}
    ICll1	= IGNORE; 		{iodine chloride}
    IBrl1	= IGNORE; 		{iodine bromide}
    CH3Il1	= IGNORE; 			{iodomethane}
    CH2I2l1	= IGNORE; 			{diiodomethane}
    CH2ClIl1	= IGNORE; 		{chloroiodomethane}
    C3H7Il1	= IGNORE; 			{2-iodopropane} 
    CH2BrIl1    = IGNORE;               {bromoiodomethane}
    CHBr2Il1    = IGNORE;               {dibromoiodomethane}
    C2H5Il1     = IGNORE;               {iodoethane}
    HIO3l1      = IGNORE;                       {}

    CH3OHl1	= 			IGNORE;
    C2H5OHl1	= 			IGNORE;
    NHSl1	= 			IGNORE;
    ClONOl1	= 			IGNORE;
    ClO3l1	= 			IGNORE;
    Cl2O3l1	= 			IGNORE;
    Br2Ol1	= 			IGNORE;
    INOl1	= 			IGNORE;
    I2O4l1	= 			IGNORE;
    I2O5l1	= 			IGNORE;
    RCll1	= 			IGNORE;
    RBrl1	= 			IGNORE;
    XORl1	= 			IGNORE;
    H2l1	= 			IGNORE;
    SORl1	= 			IGNORE;
    SPANl1	= 			IGNORE;

    OHl1	= IGNORE;		{hydroxyl radical}   
    HO2l1	= IGNORE;                       {perhydroxyl radical} 
    CH3OOl1	= IGNORE;		{ CH3O2 }
    HIO2l1      = IGNORE;
    IOl1	= IGNORE;		{iodine monoxide radical}
    O2l1	= IGNORE;      
    CO2l1	= IGNORE;
    Cll1	= IGNORE;
    Brl1	= IGNORE;
    HCOOHl1     = IGNORE;
    CH3OOHl1    = IGNORE;

    DOMl1       = IGNORE;   {unspecified dissolved organic matter}


  {ions, "el"= NEGATIVE charge}

    Hpl1         =   IGNORE;  
    OHml1        =   IGNORE;  
    O2ml1        =   IGNORE; 
    NH4pl1       =   IGNORE;
    NO3ml1       =   IGNORE; 
    NO4ml1       =   IGNORE; 
    HCO3ml1      =   IGNORE; 
    HCOOml1	 =   IGNORE;
    Clml1        =   IGNORE;
    Cl2ml1       =   IGNORE; 
    ClOml1       =   IGNORE;    
    Brml1        =   IGNORE; 
    Br2ml1	 =   IGNORE;
    BrOml1       =   IGNORE; 
    BrCl2ml1     =   IGNORE; 
    Br2Clml1     =   IGNORE; 
    HSO3ml1      =   IGNORE;
    SO32ml1      =   IGNORE; 
    HSO4ml1      =  IGNORE; 
    SO42ml1      =   IGNORE;
    NO2ml1       =   IGNORE;
    CH2OHSO3ml1  =   IGNORE;
    CO3ml1       =   IGNORE;
    SO4ml1       =   IGNORE;
    ClOHml1      =   IGNORE;   
    BrOHml1      =   IGNORE;
    CH3SO2ml1    =   IGNORE; {MSI-, CH3S(O)O-}
    CH3SO3ml1    =   IGNORE; {MS-,  CH3S(OO)O-}
    HSO5ml1      =   IGNORE;
    SO3ml1       =   IGNORE;
    SO5ml1       =   IGNORE;
    Iml1         =   IGNORE;
    IO2ml1       =   IGNORE;     
    IO3ml1       =   IGNORE;
    ICl2ml1      =   IGNORE;     
    IBr2ml1      =   IGNORE;
    IClBrml1     =   IGNORE;

 {_l2: large aerosol}
    NOl2	= IGNORE;
    NO2l2	= IGNORE;
    HNO3l2	= IGNORE;
    NH3l2	= IGNORE;
    SO2l2	= IGNORE;
    SO4l2	= IGNORE;
    O3l2	= IGNORE;

    ACO2l2	= IGNORE;
    ACTAl2	= IGNORE;             		{ CH3COOH }
    HCHOl2	= IGNORE;
    ALD2l2	= 			IGNORE;  
    H2O2l2	= IGNORE;
    ROOHl2	= 			IGNORE;   
    HONOl2	= IGNORE;

    HCll2	= IGNORE;

    N2O5l2	= IGNORE;
    HNO4l2	= IGNORE;
    NO3l2	= IGNORE;
    DMSl2	= IGNORE;
    DMSOl2	= IGNORE; 		{DMSO: CH3SOCH3}
    CH3SO2l2	= IGNORE;		{}
    CH3SO3l2	= IGNORE; 		{}
    CH3SO3Hl2	= IGNORE; 		{}
    DMOOl2	= IGNORE;		{CH3SCH2OO}
    CH3Sl2	= IGNORE;		{}
    CH3SOl2	= IGNORE;		{}
    CH3SO2Hl2	= IGNORE;		{MSIA}
    DMSO2l2	= IGNORE;		{}

    COl2	= IGNORE; 			{carbon monoxide}
    HOCll2	= IGNORE; 			{hydrochlorous acid}
    ClNO2l2	= IGNORE; 			{chlorine nitrite}
    ClNO3l2	= IGNORE; 			{chlorine nitrate}
    Cl2l2	= IGNORE; 				{molecular chlorine}
    HBrl2	= IGNORE; 			{hydrogen bromide}
    HOBrl2	= IGNORE; 			{hydrobromous acid}
    BrNO2l2	= IGNORE; 		{bromine nitrite}
    BrNO3l2	= IGNORE; 			{bromine nitrate}
    Br2l2	= IGNORE; 			{molecular bromine}
    BrCll2	= IGNORE; 			{bromine chloride}
    HIl2	= IGNORE; 			{hydrogen iodide}
    HOIl2	= IGNORE;			{hypoiodous acid}
    I2O2l2	= IGNORE;			{}
    INO2l2	= IGNORE; 			{iodine nitrite}
    INO3l2	= IGNORE; 			{iodine nitrate}
    I2l2	= IGNORE; 			{molecular iodine}
    ICll2	= IGNORE; 			{iodine chloride}
    IBrl2	= IGNORE; 			{iodine bromide}
    CH3Il2	= IGNORE; 			{iodomethane}
    CH2I2l2	= IGNORE; 			{diiodomethane}
    CH2ClIl2	= IGNORE; 		{chloroiodomethane}
    C3H7Il2	= IGNORE; 		{2-iodopropane} 
    CH2BrIl2    = IGNORE;               {bromoiodomethane} 
    CHBr2Il2    = IGNORE;               {dibromoiodomethane}
    C2H5Il2     = IGNORE;               {iodoethane}
    HIO3l2      = IGNORE;                       {}

    CH3OHl2	= 			IGNORE;
    C2H5OHl2	= 			IGNORE;
    NHSl2	= 			IGNORE;
    ClONOl2	= 			IGNORE;
    ClO3l2	= 			IGNORE;
    Cl2O3l2	= 			IGNORE;
    Br2Ol2	= 			IGNORE;
    INOl2	= 			IGNORE;
    I2O4l2	= 			IGNORE;
    I2O5l2	= 			IGNORE;
    RCll2	= 			IGNORE;
    RBrl2	= 			IGNORE;
    XORl2	= 			IGNORE;
    H2l2	= 			IGNORE;
    SORl2	= 			IGNORE;
    SPANl2	= 			IGNORE;


    OHl2	= IGNORE;
    HO2l2	= IGNORE;
    CH3OOl2	= IGNORE;
    HIO2l2      = IGNORE;
    IOl2	= IGNORE;
    O2l2	= IGNORE;
    CO2l2	= IGNORE;
    Cll2	= IGNORE;
    Brl2	= IGNORE;
    HCOOHl2     = IGNORE;
    CH3OOHl2    = IGNORE;

    DOMl2       = IGNORE;   {unspecified dissolved organic matter}

  {ions, "el"= NEGATIVE charge}

    Hpl2         =   IGNORE;
    OHml2        =   IGNORE;
    O2ml2        =   IGNORE;
    NH4pl2       =   IGNORE;
    NO3ml2       =   IGNORE;
    NO4ml2       =   IGNORE;
    HCO3ml2      =   IGNORE;
    HCOOml2	 =   IGNORE;
    Clml2        =   IGNORE;
    Cl2ml2       =   IGNORE;
    ClOml2       =   IGNORE;
    Brml2        =   IGNORE;
    Br2ml2	 =   IGNORE;
    BrOml2       =   IGNORE;
    BrCl2ml2     =   IGNORE;
    Br2Clml2     =   IGNORE;
    HSO3ml2      =   IGNORE;
    SO32ml2      =   IGNORE;
    HSO4ml2      =   IGNORE;
    SO42ml2      =   IGNORE;
    NO2ml2       =   IGNORE;
    CH2OHSO3ml2  =   IGNORE;
    CO3ml2       =   IGNORE;
    SO4ml2       =   IGNORE;
    ClOHml2      =   IGNORE;   
    BrOHml2      =   IGNORE;
    CH3SO2ml2    =   IGNORE;
    CH3SO3ml2    =   IGNORE;
    HSO5ml2      =   IGNORE;
    SO3ml2       =   IGNORE;
    SO5ml2       =   IGNORE;
    Iml2         =   IGNORE;
    IO2ml2       =   IGNORE;     
    IO3ml2       =   IGNORE;
    ICl2ml2      =   IGNORE;     
    IBr2ml2      =   IGNORE;
    IClBrml2     =   IGNORE;


{liquid phase}
 {_l3: small droplets}
    NOl3	= IGNORE;			{nitric oxide}
    NO2l3	= IGNORE;			{nitrogen dioxide}
    HNO3l3	= IGNORE;			{nitric acid}
    NH3l3	= IGNORE;			{N + 3Hammonia}
    SO2l3	= IGNORE;			{sulfur dioxide}
    SO4l3	= IGNORE;  
    O3l3	= IGNORE;			{ozone}

    ACO2l3	= IGNORE; 			{ HCOOH }
    ACTAl3	= IGNORE;             		{ CH3COOH }
    HCHOl3	= IGNORE;			{formaldehyde}
    ALD2l3	= 			IGNORE;
    H2O2l3	= IGNORE;		{hydrogen peroxide}
    ROOHl3	= 			IGNORE;   
    HONOl3	= IGNORE;		{nitrous acid}

    HCll3	= IGNORE;			{hydrogen chloride}

    N2O5l3	= IGNORE;		{dinitrogen pentoxide}    
    HNO4l3	= IGNORE;			{pernitric acid}
    NO3l3	= IGNORE;		{nitrogen trioxide}                            
    DMSl3	= IGNORE;			{ CH3-S-CH3 }   
    DMSOl3	= IGNORE; 		{DMSO: CH3SOCH3}
    CH3SO2l3	= IGNORE; 		{}
    CH3SO3l3	= IGNORE; 		{}
    CH3SO3Hl3	= IGNORE; 		{}
    DMOOl3	= IGNORE;		{CH3SCH2OO}
    CH3Sl3	= IGNORE;		{}
    CH3SOl3	= IGNORE;		{}
    CH3SO2Hl3	= IGNORE;		{MSIA}
    DMSO2l3	= IGNORE;		{}

    COl3	= IGNORE; 			{carbon monoxide}
    HOCll3	= IGNORE;			{hydrochlorous acid}
    ClNO2l3	= IGNORE;			{chlorine nitrite}
    ClNO3l3	= IGNORE;			{chlorine nitrate}
    Cl2l3	= IGNORE;				{molecular chlorine}
    HBrl3	= IGNORE; 			{hydrogen bromide}
    HOBrl3	= IGNORE;			{hydrobromous acid}
    BrNO2l3	= IGNORE; 			{bromine nitrite}
    BrNO3l3	= IGNORE; 			{bromine nitrate}
    Br2l3	= IGNORE; 				{molecular bromine}
    BrCll3	= IGNORE; 			{bromine chloride}
    HIl3	= IGNORE; 			{hydrogen iodide}
    HOIl3	= IGNORE; 			{hypoiodous acid}
    I2O2l3	= IGNORE; 			{}
    INO2l3	= IGNORE; 			{iodine nitrite}
    INO3l3	= IGNORE; 		{iodine nitrate}
    I2l3	= IGNORE; 				{molecular iodine}
    ICll3	= IGNORE; 		{iodine chloride}
    IBrl3	= IGNORE; 		{iodine bromide}
    CH3Il3	= IGNORE; 			{iodomethane}
    CH2I2l3	= IGNORE; 			{diiodomethane}
    CH2ClIl3	= IGNORE; 		{chloroiodomethane}
    C3H7Il3	= IGNORE; 			{2-iodopropane} 
    CH2BrIl3    = IGNORE;               {bromoiodomethane} 
    CHBr2Il3    = IGNORE;               {dibromoiodomethane}
    C2H5Il3     = IGNORE;               {iodoethane}
    HIO3l3      = IGNORE;                       {}

    CH3OHl3	= 			IGNORE;
    C2H5OHl3	= 			IGNORE;
    NHSl3	= 			IGNORE;
    ClONOl3	= 			IGNORE;
    ClO3l3	= 			IGNORE;
    Cl2O3l3	= 			IGNORE;
    Br2Ol3	= 			IGNORE;
    INOl3	= 			IGNORE;
    I2O4l3	= 			IGNORE;
    I2O5l3	= 			IGNORE;
    RCll3	= 			IGNORE;
    RBrl3	= 			IGNORE;
    XORl3	= 			IGNORE;
    H2l3	= 			IGNORE;
    SORl3	= 			IGNORE;
    SPANl3	= 			IGNORE;



    OHl3	= IGNORE;		{hydroxyl radical}   
    HO2l3	= IGNORE;                       {perhydroxyl radical} 
    CH3OOl3	= IGNORE;		{ CH3O2 }
    HIO2l3      = IGNORE;
    IOl3	= IGNORE;		{iodine monoxide radical}
    O2l3	= IGNORE;      
    CO2l3	= IGNORE;
    Cll3	= IGNORE;
    Brl3	= IGNORE;
    HCOOHl3     = IGNORE;
    CH3OOHl3    = IGNORE;

    DOMl3       = IGNORE;   {unspecified dissolved organic matter}


  {ions, "el"= NEGATIVE charge}

    Hpl3         =   IGNORE;  
    OHml3        =   IGNORE;  
    O2ml3        =   IGNORE; 
    NH4pl3       =   IGNORE;
    NO3ml3       =   IGNORE; 
    NO4ml3       =   IGNORE; 
    HCO3ml3      =   IGNORE; 
    HCOOml3	 =   IGNORE;
    Clml3        =   IGNORE;
    Cl2ml3       =   IGNORE; 
    ClOml3       =   IGNORE;    
    Brml3        =   IGNORE; 
    Br2ml3	 =   IGNORE;
    BrOml3       =   IGNORE; 
    BrCl2ml3     =   IGNORE; 
    Br2Clml3     =   IGNORE; 
    HSO3ml3      =   IGNORE;
    SO32ml3      =   IGNORE; 
    HSO4ml3      =  IGNORE; 
    SO42ml3      =   IGNORE;
    NO2ml3       =   IGNORE;
    CH2OHSO3ml3  =   IGNORE;
    CO3ml3       =   IGNORE;
    SO4ml3       =   IGNORE;
    ClOHml3      =   IGNORE;   
    BrOHml3      =   IGNORE;
    CH3SO2ml3    =   IGNORE;
    CH3SO3ml3    =   IGNORE;
    HSO5ml3      =   IGNORE;
    SO3ml3       =   IGNORE;
    SO5ml3       =   IGNORE;
    Iml3         =   IGNORE;
    IO2ml3       =   IGNORE;     
    IO3ml3       =   IGNORE;
    ICl2ml3      =   IGNORE;     
    IBr2ml3      =   IGNORE;
    IClBrml3     =   IGNORE;



{liquid phase}
 {_l4: large droplets}
    NOl4	= IGNORE;			{nitric oxide}
    NO2l4	= IGNORE;			{nitrogen dioxide}
    HNO3l4	= IGNORE;			{nitric acid}
    NH3l4	= IGNORE;			{N + 3Hammonia}
    SO2l4	= IGNORE;			{sulfur dioxide}
    SO4l4	= IGNORE;  
    O3l4	= IGNORE;			{ozone}

    ACO2l4	= IGNORE; 			{ HCOOH }
    ACTAl4	= IGNORE;             		{ CH3COOH }
    HCHOl4	= IGNORE;			{formaldehyde}
    ALD2l4	= 			IGNORE;  
    H2O2l4	= IGNORE;		{hydrogen peroxide}
    ROOHl4	= 			IGNORE;   
    HONOl4	= IGNORE;		{nitrous acid}

    HCll4	= IGNORE;			{hydrogen chloride}

    N2O5l4	= IGNORE;		{dinitrogen pentoxide}    
    HNO4l4	= IGNORE;			{pernitric acid}
    NO3l4	= IGNORE;		{nitrogen trioxide}                            
    DMSl4	= IGNORE;			{ CH3-S-CH3 }   
    DMSOl4	= IGNORE; 		{DMSO: CH3SOCH3}
    CH3SO2l4	= IGNORE; 		{}
    CH3SO3l4	= IGNORE; 		{}
    CH3SO3Hl4	= IGNORE; 		{}
    DMOOl4	= IGNORE;		{CH3SCH2OO}
    CH3Sl4	= IGNORE;		{}
    CH3SOl4	= IGNORE;		{}
    CH3SO2Hl4	= IGNORE;		{MSIA}
    DMSO2l4	= IGNORE;		{}

    COl4	= IGNORE; 			{carbon monoxide}
    HOCll4	= IGNORE;			{hydrochlorous acid}
    ClNO2l4	= IGNORE;			{chlorine nitrite}
    ClNO3l4	= IGNORE;			{chlorine nitrate}
    Cl2l4	= IGNORE;				{molecular chlorine}
    HBrl4	= IGNORE; 			{hydrogen bromide}
    HOBrl4	= IGNORE;			{hydrobromous acid}
    BrNO2l4	= IGNORE; 			{bromine nitrite}
    BrNO3l4	= IGNORE; 			{bromine nitrate}
    Br2l4	= IGNORE; 				{molecular bromine}
    BrCll4	= IGNORE; 			{bromine chloride}
    HIl4	= IGNORE; 			{hydrogen iodide}
    HOIl4	= IGNORE; 			{hypoiodous acid}
    I2O2l4	= IGNORE; 			{}
    INO2l4	= IGNORE; 			{iodine nitrite}
    INO3l4	= IGNORE; 		{iodine nitrate}
    I2l4	= IGNORE; 				{molecular iodine}
    ICll4	= IGNORE; 		{iodine chloride}
    IBrl4	= IGNORE; 		{iodine bromide}
    CH3Il4	= IGNORE; 			{iodomethane}
    CH2I2l4	= IGNORE; 			{diiodomethane}
    CH2ClIl4	= IGNORE; 		{chloroiodomethane}
    C3H7Il4	= IGNORE; 			{2-iodopropane} 
    CH2BrIl4    = IGNORE;               {bromoiodomethane} 
    CHBr2Il4    = IGNORE;               {dibromoiodomethane}
    C2H5Il4     = IGNORE;               {iodoethane}
    HIO3l4      = IGNORE;                       {}

    CH3OHl4	= 			IGNORE;
    C2H5OHl4	= 			IGNORE;
    NHSl4	= 			IGNORE;
    ClONOl4	= 			IGNORE;
    ClO3l4	= 			IGNORE;
    Cl2O3l4	= 			IGNORE;
    Br2Ol4	= 			IGNORE;
    INOl4	= 			IGNORE;
    I2O4l4	= 			IGNORE;
    I2O5l4	= 			IGNORE;
    RCll4	= 			IGNORE;
    RBrl4	= 			IGNORE;
    XORl4	= 			IGNORE;
    H2l4	= 			IGNORE;
    SORl4	= 			IGNORE;
    SPANl4	= 			IGNORE;


    OHl4	= IGNORE;		{hydroxyl radical}   
    HO2l4	= IGNORE;                       {perhydroxyl radical} 
    CH3OOl4	= IGNORE;		{ CH3O2 }
    HIO2l4      = IGNORE;
    IOl4	= IGNORE;		{iodine monoxide radical}
    O2l4	= IGNORE;      
    CO2l4	= IGNORE;
    Cll4	= IGNORE;
    Brl4	= IGNORE;
    HCOOHl4     = IGNORE;
    CH3OOHl4    = IGNORE;

    DOMl4       = IGNORE;   {unspecified dissolved organic matter}


  {ions, "el"= NEGATIVE charge}

    Hpl4         =   IGNORE;  
    OHml4        =   IGNORE;  
    O2ml4        =   IGNORE; 
    NH4pl4       =   IGNORE;
    NO3ml4       =   IGNORE; 
    NO4ml4       =   IGNORE; 
    HCO3ml4      =   IGNORE; 
    HCOOml4	 =   IGNORE;
    Clml4        =   IGNORE;
    Cl2ml4       =   IGNORE; 
    ClOml4       =   IGNORE;    
    Brml4        =   IGNORE; 
    Br2ml4	 =   IGNORE;
    BrOml4       =   IGNORE; 
    BrCl2ml4     =   IGNORE; 
    Br2Clml4     =   IGNORE; 
    HSO3ml4      =   IGNORE;
    SO32ml4      =   IGNORE; 
    HSO4ml4      =  IGNORE; 
    SO42ml4      =   IGNORE;
    NO2ml4       =   IGNORE;
    CH2OHSO3ml4  =   IGNORE;
    CO3ml4       =   IGNORE;
    SO4ml4       =   IGNORE;
    ClOHml4      =   IGNORE;   
    BrOHml4      =   IGNORE;
    CH3SO2ml4    =   IGNORE;
    CH3SO3ml4    =   IGNORE;
    HSO5ml4      =   IGNORE;
    SO3ml4       =   IGNORE;
    SO5ml4       =   IGNORE;
    Iml4         =   IGNORE;
    IO2ml4       =   IGNORE;     
    IO3ml4       =   IGNORE;
    ICl2ml4      =   IGNORE;     
    IBr2ml4      =   IGNORE;
    IClBrml4     =   IGNORE;




#DEFRAD                                                  
    OH		= O + H;			{hydroxyl radical}   
    HO2		= H + 2O;                       {perhydroxyl radical} 
    AHO2	= 3H + 3O + C; 			{ HOCH2O2 }    
    MCO3	= 2C + 3H + 3O; 		{ CH3CO3 }   
    MO2		= C + 3H + 2O;			{ CH3O2 }     
    ETO2	= 2C + 5H + 2O;			{ C2H5O2 }   
    KO2		= 			IGNORE;                        
    R3O2	= 			IGNORE;    
    RAO2	= 			IGNORE;   
    TO2		= 			IGNORE;     
    TCO3	= 4C + 3H + 4O; 		{ CHOCH=CHCO3 }  
    ZO2		= 			IGNORE;                        
    EO2		= 			IGNORE; { ethanolperoxy H2C(OH)CH2OO }
    PO2		= 			IGNORE;
    CHO2	= C + 2H + 2O;			{ CH2O2 / Criegee biradical}
    CRO2	= 2C + 4H + 2O;			{ CH3CHO2 }
    PRN1	= N + IGNORE;
    O1D		= O ;				{oxygen atomic first singlet state}
    O3P		= O ;				{oxygen atomic ground state}
    Cl		= Cl ;				{chlorine atomic ground state (2P3/2)}
    ClO		= Cl + O ;			{chlorine monoxide radical}
    OClO	= Cl + 2O ;			{symmetrical chlorine dioxide}
    Br		= Br ;				{bromine atomic ground state (2P3/2)}
    BrO		= Br + O ;			{bromine monoxide radical}
    I		= I ;				{iodine atomic ground state}
    IO		= I + O ;			{iodine monoxide radical}    
    OIO         = I + 2O  ;                     {iodine dioxide}

    ClRO2       = 			IGNORE; {generic chlorinated peroxy radical}
    BrRO2       = 			IGNORE; {generic brominated peroxy radical}
    IRO2        = 			IGNORE; {generic iodinated peroxy radical}



#DEFFIX
    O2		= 2O ;
    H2O		= 2H + O ;
    N2		= 2N ;

    H2Ol1	= 2H + O ;
    H2Ol2	= 2H + O ;
    H2Ol3	= 2H + O ;
    H2Ol4	= 2H + O ;
