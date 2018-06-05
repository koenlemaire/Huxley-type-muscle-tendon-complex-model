# Huxley-type-muscle-tendon-complex-model
Tutorial verson of Huxley type muscle tendon complex model. 

The current Malab implementation provides a tutorial version of the model described in the reference below. Running the tutorial will run a simulation of a massless muscle tendon complex (MTC) suspended from the ceiling, with a mass attached to the end. The MTC consists of a Huxley type contractile element, which is both in parallel and in series to elastic elements, which are simple quadratic springs. The simulation is somewhat robust against variation in parameter values and inputs, although not all feasible inputs will automatically lead to correct output, as the initial state should depend on model inputs. Please see documentation or ask the author for additional detail.

Reference:

Lemaire, K. K., Baan, G. C., Jaspers, R. T., & van Soest, A. J. (2016). Comparison of the validity of Hill and Huxley muscle tendon complex models using experimental data obtained from rat m. Soleus in situ. Journal of Experimental Biology, 219, 977-987. DOI: 10.1242/jeb.128280

-- main script (will run as is and produce model output, provided that all other functions are on the matlab path)

run_hux_tutorial.m

-- main function

hux_tutorial.m

All other functions are subfunctions to compute initial conditions and state derivatives. 

All files are released under the terms of the GNU General Public License,
version 3. See http://www.gnu.org/licenses/gpl.html

--------------------------------------------------------------------------

Author: Koen Lemaire (kklemaire_edu@posteo.nl)
