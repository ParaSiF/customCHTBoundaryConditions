/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file NeumannCoupledBoundary.H
 * @author B McCormick
 * @date 3 March 2026
 * @brief Header file for NeumannCoupledBoundary.C
 */

#ifndef NeumannCoupledBoundary_H
#define NeumannCoupledBoundary_H

#include "../CoupledBoundary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
   Class NeumannCoupledBoundary Declaration
\*---------------------------------------------------------------------------*/

class NeumannCoupledBoundary
:
    public CoupledBoundary
{

public:

    //- Runtime type information
    TypeName("NeumannCoupledBoundary");


    // Constructors

        //- Construct from patch and internal field
        NeumannCoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        NeumannCoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  NeumannCoupledBoundary onto a
        //  new patch
        NeumannCoupledBoundary
        (
            const NeumannCoupledBoundary&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        NeumannCoupledBoundary
        (
            const NeumannCoupledBoundary&
        );

        //- Construct as copy setting internal field reference
        NeumannCoupledBoundary
        (
            const NeumannCoupledBoundary&,
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
