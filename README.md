gbc-boxmodel
============

Amos et al. (2013; 2014) box model for global Hg cycling

Last updated: 04 September 2017

Questions or comments? Contact:
  Elsie Sunderland or Colin Thackray
  Harvard University
  
  Email: ems@seas.harvard.edu or thackray@seas.harvard.edu
  Web: http://bgc.seas.harvard.edu

## Citation for code
 1. Amos, H. M., D.J. Jacob, D.G. Streets, E.M. Sunderland (2013), Legacy impacts of all-time anthropogenic emissions on the global mercury cycle, Glob. Biogeochem. Cycle, 27(2), 410-421.
 1. Amos, H. M., et al. (2014), Global biogeochemical implications of mercury discharges from rivers and sediment burial, Environ. Sci.Technol., 48(16), 9514-9522.

## Credit

 * Co-authorship is appropriate if your paper benefits significantly from use of this model/code.
 * Citation is appropriate if use of this model/code has only a marginal impact on your work or if the work is a second generation application of the model/code.

### Submitting bugs or science updates

If you find a bug, please report it to us (thackray@seas.harvard.edu)
with the subject "Hg box model: bug report". We'll fix it, document
it, and post corrected code online.  If you'd like to submit a science
update, please contact me (thackray@seas.harvard.edu) with the subject
"Hg box model: science update". In the email, provide a quick
description of the update and a copy of the journal article associated
with the update. We will merge the update into the standard version of
the code available online.

We're excited to include science updates from the community, but as a
policy we will ONLY include work associated with a published/accepted
manuscript.  Thanks for your participation and interest!

## Software requirements

The box model is written in MATLAB version R2012b. The code will
probably run on older versions, but hasn’t been tested. If the code
does not run on your version of MATLAB, please let us know.

## Code version

This release is code version 2.0 and corresponds the updated box model
published in Amos et al. (2014, ES&T).


## Running the code
 1. Open MATLAB.
 1. Open forWeb_main.m. This module is the driver for the whole box module.
 1. Type the following in the Command Window:
     > forWeb_main

or press the green “Run” button: Press the green “Run” button

## Model output

When you run the code for the first time, a few plots will pop up and
model output will print to the Command Window. The plots and output
are intended to give the user a basic sense of what the model is
doing. If you’d like to turn off the plots, set Lplot = 0 in
forWeb_main.m. If you’d like to suppress output to the Command Window,
set Ldisp = 0 in forWeb_main.m.


## Code structure

The code is highly modular. Each module performs a distinct task and
they’re called in sequence from forWeb_main.m. Briefly, the modules
are separated into handling anthropogenic emissions, setting up rate
coefficients for solving the coupled system of ODEs, and running the
preanthropogenic and anthropogenic simulations.  

`forWeb_main.m`: This is the main driver for the box model. Users can
select different options for anthropogenic emission inventories,
discharge of Hg from rivers, or future vs. presentday simulations.

`forWeb_AnthroEmiss.m`: Sets up anthropogenic Hg emissions.

`forWeb_product_emissions.m`; Called from forWeb_AnthroEmiss.m. This
module deals with emissions from intentional use of Hg in products and
processes (Horowitz et al., 2014).

`forWeb_rate_coefficients.m`: Calculates first-order rate
coefficients. The cycling of Hg among the ocean, atmosphere, and
terrestrial ecosystems is represented as a set of coupled first-order
differential equations, which depend on first order rate coefficients.

`forWeb_makeA.m`: Assembles first-order rate coefficients into matrix
A, known as the Jacobian, which is used to solved the system of
coupled ODEs.

`forWeb_RunPreAnthro.m`: Simulates the pre-anthropogenic natural
steady-state of mercury, which is used as the initial conditions for
simulating the anthropogenic era.

`forWeb_RunAnthro.m`: Simulations the anthropogenic era, forced by
whichever anthropogenic emission inventory the user has selected in
forWeb_main.m.

`forWeb_displayOutput.m`: Prints model results to the Command Window.
