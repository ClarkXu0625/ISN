Folder Single_Unit_Sweep tests a E-I coupled unit system's bistability, 
under which parameters, including different connections and thresholds, are 
sweeped. Cross connections is not considered in this session.

'function' folder includes all the generic testing and simulation functions 
that helps modify parameter anlsysis uniformly.

'extended_range' folder expands the parameter range to finder trends in 
larger scale, but also requires significantly longer running time.

'EE_EI.m' sweeps EE and EI connection strengths to config unit bistability
property.

'EE_EI_others.m' gives different EE-EI plots under different WII connection
strengths, in order to find the trend in higher parameter dimension