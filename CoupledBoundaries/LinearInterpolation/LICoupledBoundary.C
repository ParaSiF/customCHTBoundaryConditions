/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file LICoupledBoundary.C
 * @author B McCormick
 * @date 3 March 2026
 * @brief Custom coupled boundary condition imposing a fixed value (temperature) on both#
 * sides of the interface based on weighted average of temperature
 */

#include "LICoupledBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LICoupledBoundary::
LICoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    CoupledBoundary(p, iF)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


LICoupledBoundary::
LICoupledBoundary
(
    const LICoupledBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    CoupledBoundary(ptf, p, iF, mapper)
{}


LICoupledBoundary::
LICoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
   CoupledBoundary(p, iF, dict)
{

    // define quantities for pushing/fetching with mui
    push_quantity = "temp";
    fetch_quantity = "flux";

    this->readValueEntry(dict, IOobjectOption::MUST_READ);
    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }
}


LICoupledBoundary::
LICoupledBoundary
(
    const LICoupledBoundary& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    CoupledBoundary(wtcsf, iF)
{}


LICoupledBoundary::
LICoupledBoundary
(
    const 
    LICoupledBoundary& wtcsf
)
:
    CoupledBoundary(wtcsf)
{}


void LICoupledBoundary::updateCoeffs()
{
    if (updated() || this->size() == 0)
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType(); 

    // get cell centre coordinates of boundary cells
    vectorField coords = patch().Cf();

    // k/dx for computing weight of average
    const scalarField KDelta = kappa(*this)*patch().deltaCoeffs();

    // push k/dx to neighbour, note that this assumes constant thermal conductivity and discretisation
    push(coords[0].x(), coords[0].y(), coords[0].z(), KDelta[0], iteration);

    // fetch k/dx from neighbour, note that this assumes constant thermal conductivity and discretisation
    scalar nbrKDelta = fetch(coords[0].x(), coords[0].y(), coords[0].z(), iteration);

    // set weight of average
    this->valueFraction() = nbrKDelta/(nbrKDelta+KDelta);

    // push internal temperature field at boundary to neighbour
    scalarField internal = patchInternalField();
    forAll(internal, faceI){
        push(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), internal[faceI], iteration);
    }  

    // fetch internal temperature field at boundary from neighbour
    scalarField nbrIntFld = scalarField(this->size());
    forAll(nbrIntFld, faceI){
        nbrIntFld[faceI] = fetch(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), iteration);
    }


    // set reference value to neighbour internal temperature field
    this->refValue() = nbrIntFld;
    this->refGrad() = 0;

    CoupledBoundary::updateCoeffs();
    UPstream::msgType(oldTag);  // Restore tag
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    
LICoupledBoundary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


// ************************************************************************* //
