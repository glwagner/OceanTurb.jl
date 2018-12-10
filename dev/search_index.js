var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#OceanBoundaryLayerModels.jl-1",
    "page": "Home",
    "title": "OceanBoundaryLayerModels.jl",
    "category": "section",
    "text": "OceanBoundaryLayerModels.jl implements models for the ocean\'s boundary layer."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "As simple as git clone https://github.com/glwagner/OceanBoundaryLayerModels.jl.git"
},

{
    "location": "#Authors-1",
    "page": "Home",
    "title": "Authors",
    "category": "section",
    "text": "Gregory L. Wagner"
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
    "text": "Models for the oceanic boundary layer are partial differential equations that approximate the effects of Internal and surface fluxes of heat, salinity, and momentum due to\nabsorption of incoming solar radiation;\ncooling by outgoing radiation;\nlatent and sensible heat exchange with the atmosphere;\nevaporation and precipitation;\nwind forcing;\nVertical turbulent fluxes due to\nconvection / gravitational instability;\nmechanical mixing due to wind and boundary current shear OceanMixedLayerModels.jl uses  an implementation of atmospheric and radiative forcings that is shared across all models. The models therefore differ mainly in the way they parameterize convective and mechanical mixing."
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
    "text": "The one-dimensional, horizontally-averaged boundary-layer equations for  horizontal momentum, salinity, and temperature are beqs\nu_t =   f v - G^u_z - F^u c labelxmomentum \nv_t = - f u - G^v_z - F^v c \nS_t =       - G^S_z - F^S c \nT_t =       - G^T_z - F^T c labeltemperature\neeqswhere subscripts t and z denote derivatives with respect to time  and the vertical coordinate z and f is the Coriolis parameter.  In \\eqref{xmomentum}–\\eqref{temperature}, internal and boundary forcing of a quantity phi is denoted F^phi, while vertical turbulent fluxes areG^phi = overlinew phi p\n\n\n Numerical modeling of the oceanic boundary layer\n\n Basic form\n\nModels for the oceanic boundary layer are partial differential equations of \nthe form\nmath \\beq \\phi_t = L \\phi + N(\\phi) \\c \\label{mathematicalform}faf0627972293e4cfa03da9de9e9daa31143d1d5:docs/src/numerics.md\\eeq\nwhere ``phi`` is a variable like velocity, temperature, or salinity, ``L`` is \na linear operator, and ``N`` is a nonlinear operator.\n\n## Time-stepping\n\nAn explicit forward Euler time integration scheme discretizes\n\\eqref{mathematicalform} in time with\nmath \\beq \\phi^{n+1} = \\phi^{n} + \\Delta t \\left [ L \\phi^n + N(\\phi^n) \\right ] \\, \\eeq  ```where the superscripts n and n+1 denote the solution at  time-step n and n+1, respectively."
},

{
    "location": "models/pricewellerpinkel/#",
    "page": "OceanMixedLayerModels.PriceWellerPinkel",
    "title": "OceanMixedLayerModels.PriceWellerPinkel",
    "category": "page",
    "text": ""
},

{
    "location": "models/pricewellerpinkel/#OceanMixedLayerModels.PriceWellerPinkel-1",
    "page": "OceanMixedLayerModels.PriceWellerPinkel",
    "title": "OceanMixedLayerModels.PriceWellerPinkel",
    "category": "section",
    "text": "This module implements a model for the oceanic boundary layer conceived by Price, Weller, and Pinkel (1986)."
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
    "text": ""
},

{
    "location": "man/types/#OceanBoundaryLayerModels.ForcingData-Tuple{}",
    "page": "Private types",
    "title": "OceanBoundaryLayerModels.ForcingData",
    "category": "method",
    "text": "Forcing(; tdata=[0, year], forcingfields...)\n\nConstruct a Forcing specified at time points in tdata. The forcing fields, which are arrays of data corresponding to the times in tdata and  are assumed to be measured/calculated/specified at the surface and z=0, are specified by keyword arguments. Their default values are 0tdata. Possible forcing inputs are\n\nshortwave : incoming shortwave radiation\nlongwave  : outgoing longwave radiation\nlatent    : incoming latent heat flux\nsensible  : incoming latent heat flux\nprecip    : precipitation\nevap      : evaporation\nxstress   : wind stress in the x-direction\nystress   : wind stress in the y-direction\n\n\n\n\n\n"
},

{
    "location": "man/types/#Private-types-in-module-OceanBoundaryLayerModels:-1",
    "page": "Private types",
    "title": "Private types in module OceanBoundaryLayerModels:",
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
    "text": ""
},

{
    "location": "man/functions/#Functions-exported-from-OceanBoundaryLayerModels:-1",
    "page": "Functions",
    "title": "Functions exported from OceanBoundaryLayerModels:",
    "category": "section",
    "text": "Modules = [OceanBoundaryLayerModels]\nPrivate = false\nOrder = [:function]"
},

]}
