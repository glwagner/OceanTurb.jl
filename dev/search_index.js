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
    "location": "physics/#",
    "page": "Physics",
    "title": "Physics",
    "category": "page",
    "text": "newcommandc \n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "physics/#Physics-of-the-oceanic-boundary-layer-1",
    "page": "Physics",
    "title": "Physics of the oceanic boundary layer",
    "category": "section",
    "text": "The dynamics of the ocean\'s boundary layer are dominated by forcing from the atmosphere. Atmospheric forcing includes momentum forcing by winds, salinity forcing by evaporation and precipitation, and heat forcing by latent heat fluxes, sensible heat fluxes, and incoming and outgoing radiation."
},

{
    "location": "numerics/#",
    "page": "Modeling",
    "title": "Modeling",
    "category": "page",
    "text": "newcommandc \n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "numerics/#Numerical-modeling-of-the-oceanic-boundary-layer-1",
    "page": "Modeling",
    "title": "Numerical modeling of the oceanic boundary layer",
    "category": "section",
    "text": "Models for the oceanic boundary layer are partial differential equations that approximate the effects of Atmospheric and radiative forcing;\nParameterization of convection due to surface cooling;\nParameterization of mechanical mixing by mixed layer turbulence due mainly to wind forcing of boundary-layer currents.OceanMixedLayerModels.jl uses  an implementation of atmospheric and radiative forcings that is shared across all models. The models therefore differ mainly in the way they parameterize convective and mechanical mixing."
},

{
    "location": "numerics/#Coordinate-system-1",
    "page": "Modeling",
    "title": "Coordinate system",
    "category": "section",
    "text": "We use a Cartesian coordinate system in which gravity points downwards,  toward the ground or bottom of the ocean. The vertical coordinate z  thus increases upwards. We locate the surface at z=0. This means that if the boundary layer has depth h, the bottom of the boundary layer is  located at z=-h."
},

{
    "location": "numerics/#Governing-equations-1",
    "page": "Modeling",
    "title": "Governing equations",
    "category": "section",
    "text": "The one-dimensional, horizontally-averaged boundary-layer equations for  horizontal momentum, salinity, and temperature are beqs\nu_t - f v = -G^u_z - F^u_z c \nv_t + f u = -G^v_z - F^v_z c \n      S_t = -G^S_z - F^S_z c \n      T_t = -G^T_z - F^T_z c\neeqswhere subscripts t and z denote derivatives with respect to time  and the vertical coordinate z, f is the Coriolis parameter,  G^phi = overlinew phi denotes the turbulent vertical flux of  a quantity phi, while F^phi denotes vertical fluxes due to  forcing."
},

{
    "location": "numerics/#Temperature-forcing-1",
    "page": "Modeling",
    "title": "Temperature forcing",
    "category": "section",
    "text": "We write the temperature forcing F^T asF^T = F^rlat + F^rsens + F^rlongwave \n        + F^rshortwave cin terms of the four contributions from latent heating, sensible heating,  outgoing longwave radiation, and incoming shortwave radiation.  The first three contributions are implemented as effective boundary conditions in the uppermost gridpoints of the model. Shortwave radiation, on the other hand, heats the interior of the boundary  layer.  We parameterize the effect of interior heating by shortwave radiation by dividing the shortwave spectrum into infrared (IR) and ultraviolet (UV) components and introducing a z-dependent absorption function such thatbeq\nF^rshortwave(z) = F^rshortwave_0\n    left ( alpha_rIR exp left  zd_rIR right  \n          + alpha_rUV exp left  zd_rUV right  right ) c\nlabelshortwaverad\neeqwhere F^rshortwave_0 is the incoming shortwave radiation at the  surface, which is provided as an input to the boundary layer model. In \\eqref{shortwaverad},  d_rIR and d_rUV are the  penetration scales of infrared and ultraviolet radiation, and alpha_rIR and alpha_rUV are the fractions of total shortwave radiation with infrared and ultraviolet wavelengths, respectively, such that  alpha_rIR + alpha_rUV = 1."
},

{
    "location": "numerics/#Salinity-forcing-1",
    "page": "Modeling",
    "title": "Salinity forcing",
    "category": "section",
    "text": "Salinity forcing F^S at the surface isbeq\nF^S(z=0) = S(z=0)(E-P) c\nlabelsalinityflux\neeqwhere S is salinity, E is evaporation, and P is precipitation. Equation \\eqref{salinityflux} is implemented as an effective boundary condition in the uppermost grid points of the model."
},

{
    "location": "numerics/#Momentum-forcing-1",
    "page": "Modeling",
    "title": "Momentum forcing",
    "category": "section",
    "text": "Momentum forcing at the surface is beq\nF^u(z=0) = fractau^xrho_0 c\nlabeluforcing\neeqwhere tau^x is the wind stress in the x-direction with units of  Nrm^2, and rho_0 is a reference density. A similar form is used for y-momentum. Equation \\eqref{uforcing} is implemented as an effective boundary condition."
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
    "location": "man/functions/#OceanBoundaryLayerModels.loadforcing-Tuple{Any}",
    "page": "Functions",
    "title": "OceanBoundaryLayerModels.loadforcing",
    "category": "method",
    "text": "loadforcing(filepath)\n\nInitalize a Forcing from the JLD2 file at filepath. The names of the file\'s  fields must correspond to the names of the Forcing fields.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#Functions-exported-from-OceanBoundaryLayerModels:-1",
    "page": "Functions",
    "title": "Functions exported from OceanBoundaryLayerModels:",
    "category": "section",
    "text": "Modules = [OceanBoundaryLayerModels]\nPrivate = false\nOrder = [:function]"
},

]}
