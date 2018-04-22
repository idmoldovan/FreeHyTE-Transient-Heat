# FreeHyTE-Transient-Heat
FreeHyTE - Transient Heat uses the temperature model of the hybrid-Trefftz finite element formulation for the solution of 2D transient heat
conduction problems with interior heat generation. It can be used, in a more general context, for the solution of 2D non-homogeneous 
parabolic boundary value problems with arbitrary source terms.

A generalized mid-point time stepping method is used for the time discretization. The solution in space is obtained approximating independently the temperature field in the domain of each element and the heat flux field on its essential boundaries. The trademark of the Trefftz elements is that the domain aprroximation is constructed on functions satisfying exactly the homogeneous form of the governing equations. The particular solution is obtained using a novel variant of the Dual Reciprocity Method. Due to the relevant physical information embedded in the domain bases, Trefftz elements are insensitive gross mesh distortions and recover well high field gradients. Moreover, as compared to the conventional formulations, the temperature and heat flux solutions are much more balanced in terms of quality and very large finite elements are affordable.

On the boundaries of the domain, the temperatures or heat fluxes can be enforced, or convective boundary conditions can be specified. All boundary conditions are enforced weakly. 

The implementation features an easy-to-use GUI, regular and irregular mesh generators, saving and loading of models, and Installation, User's and Developer's Manuals. The output includes plots of temperature and heat flux fields at selected instants, and data files with the results for posterior post-processing (if requested by the user).
