/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file DirichletCoupledBoundary.C
 * @author B McCormick
 * @date 3 March 2026
 * @brief Custom coupled boundary condition imposing a fixed value (temperature) condition while
 * sending heat flux data.
 */

#include "DirichletCoupledBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DirichletCoupledBoundary::
DirichletCoupledBoundary
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


DirichletCoupledBoundary::
DirichletCoupledBoundary
(
    const DirichletCoupledBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    CoupledBoundary(ptf, p, iF, mapper)
{}


DirichletCoupledBoundary::
DirichletCoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
   CoupledBoundary(p, iF, dict)
{

    // define quantities for pushing/fetching with mui
    push_quantity = "flux";
    fetch_quantity = "temp";

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


DirichletCoupledBoundary::
DirichletCoupledBoundary
(
    const DirichletCoupledBoundary& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    CoupledBoundary(wtcsf, iF)
{}


DirichletCoupledBoundary::
DirichletCoupledBoundary
(
    const 
    DirichletCoupledBoundary& wtcsf
)
:
    CoupledBoundary(wtcsf)
{}


void DirichletCoupledBoundary::updateCoeffs()
{
    if (updated() || this->size() == 0)
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType(); 
    
    // Thermal conductivity field at boundary
    const scalarField kappaTp = kappa(*this);

    // Calculate heatflux normal to boundary
    scalarField Q = -kappaTp*patch().magSf()*snGrad();
    
    // get cell centre coordinates of boundary cells
    vectorField coords = patch().Cf();
    
    // push heat flux field at boundary to neighbour
    forAll(Q, faceI){
        push(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), Q[faceI], iteration);
    }
    
    // fetch temperature at boundary from neighbour
    scalarField nbrTempFld = scalarField(this->size());
    forAll(nbrTempFld, faceI){
        nbrTempFld[faceI] = fetch(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), iteration);
    }


    // set reference value of mixed boundary conditions
    this->refGrad() = 0;
    this->refValue() = nbrTempFld;
    this->valueFraction() = 1.0;

    CoupledBoundary::updateCoeffs();
    UPstream::msgType(oldTag);  // Restore tag
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    
DirichletCoupledBoundary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


// ************************************************************************* //
