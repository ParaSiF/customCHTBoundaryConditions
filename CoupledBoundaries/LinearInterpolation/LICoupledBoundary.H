/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file LICoupledBoundary.H
 * @author B McCormick
 * @date 3 March 2026
 * @brief Header file for LICoupledBoundary.C
 */

#ifndef LICoupledBoundary_H
#define LICoupledBoundary_H

#include "../CoupledBoundary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
   Class LICoupledBoundary Declaration
\*---------------------------------------------------------------------------*/

class LICoupledBoundary
:
    public CoupledBoundary
{

public:

    //- Runtime type information
    TypeName("LinearInterpolationCoupledBoundary");


    // Constructors

        //- Construct from patch and internal field
        LICoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        LICoupledBoundary
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  LICoupledBoundary onto a
        //  new patch
        LICoupledBoundary
        (
            const LICoupledBoundary&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        LICoupledBoundary
        (
            const LICoupledBoundary&
        );

        //- Construct as copy setting internal field reference
        LICoupledBoundary
        (
            const LICoupledBoundary&,
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
