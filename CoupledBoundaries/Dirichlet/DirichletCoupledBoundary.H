/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file DirichletCoupledBoundary.H
 * @author B McCormick
 * @date 3 March 2026
 * @brief Header file for DirichletCoupledBoundary.C
 */

#ifndef DirichletCoupledBoundary_H
#define DirichletCoupledBoundary_H

#include "../CoupledBoundary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
   Class DirichletCoupledBoundary Declaration
\*---------------------------------------------------------------------------*/

class DirichletCoupledBoundary
:
    public CoupledBoundary
{

public:

    //- Runtime type information
    TypeName("DirichletCoupledBoundary");


    // Constructors

        //- Construct from patch and internal field
        DirichletCoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        DirichletCoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  DirichletCoupledBoundary onto a
        //  new patch
        DirichletCoupledBoundary
        (
            const DirichletCoupledBoundary&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        DirichletCoupledBoundary
        (
            const DirichletCoupledBoundary&
        );

        //- Construct as copy setting internal field reference
        DirichletCoupledBoundary
        (
            const DirichletCoupledBoundary&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
