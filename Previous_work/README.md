Unfinished ISN project by Connor Zawacki (connorzawacki@brandeis.edu)

# Brief Abstract/Introduction
An attractor network is a model for a system of neurons whose activity 
tends towards one of several discrete steady stable states, called 
attractor states. A network of neurons with a large number of attractor 
states with comparatively few neurons could provide a plausible framework 
for information processing. In the past, attempts to computationally 
replicate such a network have yielded inadequate results, characterized by 
problems such as unbiologically high firing rates or requiring the fine 
tuning of many parameters. There is strong experimental evidence in 
cortical and hippocampal circuits for system dynamics to follow an 
inhibition stabilized regime, characterized by strong recurrent excitation 
which is stabilized by stronger recurrent inhibition. Use of firing rate 
models in the inhibition stabilized regime yields more promising results, 
as multiple attractor states can be stabilized with arbitrary low firing 
rates. By simulating large variations of connection strengths while fixing 
other parameters, we can test both the number of possible steady states for 
a given parameter set, and potentially define a relationship between 
connection strengths required to produce multistability. We successfully 
produced several networks capable of a combinatorial explosion of attractor 
states. Moreover, this network continued to produce this large number of 
attractor states for a variety of connection strengths between units. There 
are several plausible steps forward from here, the goal of course being to 
define under exactly what conditions this combinatorial explosion of 
attractor states are stable.

# Code
Naming convention: 
ISN_net_excitFiringRate_InhibFiringRate_crossConnection_#ofPairs

Note concerning HPCC results: I specifically did not include code I have 
written for the Brandeis HPCC. The reason for this is twofold. Firstly I 
was unable to locate an error or two in what I have written for it causing 
some figures to be inaccurate. Secondly some of the code isn’t widely 
applicable (i.e. is specific to the batch scripts I have written and the 
pathnames on my machine)


single_unit_sweep.m
A single E-I pair param sweep for excitatory and Inhibitory threshold, produces a grid
plot of param combinations resulting in bistability for the single E-I pairing

single_unit_Bistable_QI.m
A single E-I pair with a nonlinear inhibitory firing rate equation exhibiting bistability

ISN_net_QE_linI_EIX_20.m
20 E-I pairs with no E-I cross connection and nonlinear firing rate equations for
excitatory units. With current parameter set anywhere between 2 and 19 of the 20 pairs can be “on” at once, and the rest are suppressed.
 
ISN_net_QE_linI_40_EIX.m
Originally written by Professor Miller. 40 total E-I pairs, 10 of which are active. Nonlinear
Firing rate equation for excitatory units, and linear for inhibitory. E-I cross connection

ISN_net_linE_QI_IEX_10.m
10 E-I pairs with linear excitatory firing rate equations and nonlinear inhibitory firing rate
equations. I-E cross connection. Anywhere from 3-5 units can be “on” at once. Worth exploring stability more, as this particular regimen seem a bit more “fragile” than the others.

ISN_net_linE_QI_0x_2.m
Originally written by Professor Miller. 2 E-I pairs, 1 of which is active. Nonlinear firing
rate equation for inhibitory units, and linear for excitatory. No cross connection, essentially just exhibiting bistability in a single unit.

ISN_net_lin_IEX_20.m
20 total E-I pairs with only linear firing rate equations. I-E cross connections. With
current parameter set anywhere between 1-18 E-I pairs can be on at once, the rest being suppressed.

ISN_net_lin_IEX_10p.m
Originally written by Professor Miller 10 E-I pairs, 8 of which are active. Linear firing rate
equations and i-e cross connection.

ISN_net_lin_IEX_10.m
Very close to ISN_net_lin_IEX_10p, just exploring slightly different parameters and
number of active E-I pairnigs.