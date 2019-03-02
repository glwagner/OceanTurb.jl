var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#OceanTurb.jl-1",
    "page": "Home",
    "title": "OceanTurb.jl",
    "category": "section",
    "text": "OceanTurb.jl implements models for the ocean\'s turbulent surface boundary layer."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "As simple as git clone https://github.com/glwagner/OceanTurb.jl.git"
},

{
    "location": "#Authors-1",
    "page": "Home",
    "title": "Authors",
    "category": "section",
    "text": "Gregory L. Wagner."
},

{
    "location": "basics/#",
    "page": "Basics",
    "title": "Basics",
    "category": "page",
    "text": "newcommandc \nnewcommandp \nnewcommanddpartial\n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "basics/#Physics-and-modeling-of-the-oceanic-boundary-layer-1",
    "page": "Basics",
    "title": "Physics and modeling of the oceanic boundary layer",
    "category": "section",
    "text": "Models for the oceanic boundary layer are partial differential equations that approximate the effects of internal and surface fluxes, and\nturbulent vertical fluxeson the evolution of temperature, salinity, and momentum in the near-surface ocean.Internal and surface fluxes of heat are due toabsorption of incoming shortwave solar radiation;\ncooling by outgoing longwave radiation;\nlatent and sensible heat exchange with the atmosphere.Surface fluxes of salinity occur due to evaporation and precipitation,  while momentum fluxes are associated with atmospheric winds.Vertical turbulent fluxes are typically associated withgravitational instability and convection, and\nmechanical turbulent mixing associated with currents and wind forcing.OceanMixedLayerModels.jl uses  an implementation of atmospheric and radiative forcings that is shared across all models. The models therefore differ mainly in the way they parameterize convective and mechanical mixing."
},

{
    "location": "basics/#Coordinate-system-1",
    "page": "Basics",
    "title": "Coordinate system",
    "category": "section",
    "text": "We use a Cartesian coordinate system in which gravity points downwards,  toward the ground or bottom of the ocean. The vertical coordinate z  thus increases upwards. We locate the surface at z=0. This means that if the boundary layer has depth h, the bottom of the boundary layer is  located at z=-h."
},

{
    "location": "basics/#Governing-equations-1",
    "page": "Basics",
    "title": "Governing equations",
    "category": "section",
    "text": "The one-dimensional, horizontally-averaged boundary-layer equations for  horizontal momentum U and V, salinity S, and  temperature T are beqs\nU_t =   f V - d_z overlinew u - F^u c labelxmomentum \nV_t = - f U - d_z overlinew v - F^v c \nS_t =       - d_z overlinew s - F^S c \nT_t =       - d_z overlinew theta - F^T c labeltemperature\neeqswhere subscripts t and z denote derivatives with respect to time  and the vertical coordinate z and f is the Coriolis parameter.  The lowercase variables u, v, s, and theta refer to the  three-dimensional perturbations from horizontal velocity, salinity, and  temperature, respectively.  In \\eqref{xmomentum}–\\eqref{temperature}, internal and boundary forcing of a quantity Phi is denoted F^Phi, while vertical turbulent fluxes areG^Phi = overlinew phi p"
},

{
    "location": "basics/#Numerics-1",
    "page": "Basics",
    "title": "Numerics",
    "category": "section",
    "text": "One of the surprising aspects of boundary layer parameterization is that  spatial discretization and time-stepping considerations are not less important or even distinct from physical aspects of approximation and modeling."
},

{
    "location": "basics/#Spatial-discretization-1",
    "page": "Basics",
    "title": "Spatial discretization",
    "category": "section",
    "text": "OceanBoundaryLayerModels.jl uses a finite volume method to discretize the oceanic boundary layer in z."
},

{
    "location": "basics/#Time-stepping-1",
    "page": "Basics",
    "title": "Time-stepping",
    "category": "section",
    "text": "Models for the oceanic boundary layer are partial differential equations of  the formbeq\nPhi_t = L Phi + N(Phi) c\nlabelmathematicalform\neeqwhere Phi is a variable like velocity, temperature, or salinity, L is  a linear operator, and N is a nonlinear operator. The linear and nonlinear parts of the boundary layer model are delineated to permit hybrid  implicit/explicit time-stepping schemes."
},

{
    "location": "basics/#Forward-Euler-method-1",
    "page": "Basics",
    "title": "Forward-Euler method",
    "category": "section",
    "text": "An explicit forward Euler time integration scheme discretizes \\eqref{mathematicalform} in time withbeq\nPhi^n+1 = Phi^n + Delta t left  L Phi^n + N(Phi^n) right  \neeq where the superscripts n and n+1 denote the solution at  time-step n and n+1, respectively."
},

{
    "location": "models/pacanowskiphilander/#",
    "page": "The Pacanowski-Philander (1981) model",
    "title": "The Pacanowski-Philander (1981) model",
    "category": "page",
    "text": "newcommandc \nnewcommandp \nnewcommanddpartial\n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "models/pacanowskiphilander/#The-Pacanowski-Philander-(1981)-model-1",
    "page": "The Pacanowski-Philander (1981) model",
    "title": "The Pacanowski-Philander (1981) model",
    "category": "section",
    "text": "[Pacanowski and Philander (1981)] propose a simple one-dimensional model for equatorial boundary layers dominated by mechanical turbulent mixing. In their model, horizontal velocity, temperature, and salinity are governed bybeqs\nU_t =   f V - d_z overlinew u - F^u c \nV_t = - f U - d_z overlinew v - F^v c \nS_t =       - d_z overlinew s - F^S c \nT_t =       - d_z overlinew theta - F^T c \neeqswhere uppercase variables are resolved, mean quantities,  and lowercase variables are unresolved perturbations. A key variable in the Pacaonwski and Philander formulation is the local Richardson number, defined byRi = - fracg rho_zrho_0 left ( U_z^2 + V_z^2 right ) cwhere g is the gravitational constant, rho_0 is a reference density, and rho is density. For simplicity, we use a linear equation of state between temperature, salinity, and density,rho(z t) = rho_0 left  \n  alpha left ( T-T_0 right ) + beta left ( S-S_0 right ) right  cwhere typical values for the temperature and salinity expansion coefficients  are alpha = 2 times 10^-4 rK m^3  kg and beta = 1."
},

{
    "location": "models/pacanowskiphilander/#Eddy-diffusivities-for-momentum,-temperature,-and-salinity-1",
    "page": "The Pacanowski-Philander (1981) model",
    "title": "Eddy diffusivities for momentum, temperature, and salinity",
    "category": "section",
    "text": "The core of the [PP81] model is to parameterize vertical turbulent fluxes  with an eddy diffusivity/eddy viscosity, such that overlinew phi = kappa_Phi d_z Phi cwhere kappa_Phi is the eddy diffusivity or viscosity of the quantity  Phi.The eddy viscosity for both U and V iskappa_U = nu_0 + fracnu_1left ( 1 + c Ri right )^n cwhile the eddy diffusivity for T and S arekappa_T = kappa_0 + frac kappa_1  left ( 1 + c Ri right )^n+1 pThis parameterization implies an Ri-dependent turbulent Prandtl number ofPr = frackappa_Ukappa_T approx left ( 1 + c Ri right ) cwhen Ri is small. Typical values for the parameters (see [CV12]) arenu_0 = 10^-4 rcms,\nnu_1 = 10^-2 rcms,\nkappa_0 = 10^-5 rcms,\nkappa_1 = 10^-2 rcms,\nc = 5, and\nn=2.[Pacanowski and Philander (1981)]: https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2) [PP81]: https://journals.ametsoc.org/doi/abs/10.1175/1520-0485(1981)011%3C1443:POVMIN%3E2.0.CO;2) [CV12]: https://books.google.com/books?id=AAfoCAAAQBAJ"
},

{
    "location": "man/types/#",
    "page": "Private types",
    "title": "Private types",
    "category": "page",
    "text": ""
},

{
    "location": "man/types/#Private-types-1",
    "page": "Private types",
    "title": "Private types",
    "category": "section",
    "text": "Modules = [OceanBoundaryLayerModels]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "man/functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [OceanBoundaryLayerModels]\nPrivate = false\nOrder = [:function]"
},

]}
