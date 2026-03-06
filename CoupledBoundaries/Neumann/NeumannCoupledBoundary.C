/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
/**
 * @file NeumannCoupledBoundary.C
 * @author B McCormick
 * @date 3 March 2026
 * @brief Custom coupled boundary condition imposing a fixed gradient (flux) condition while
 * sending temperature data.
 */

#include "NeumannCoupledBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NeumannCoupledBoundary::
NeumannCoupledBoundary
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


NeumannCoupledBoundary::
NeumannCoupledBoundary
(
    const NeumannCoupledBoundary& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    CoupledBoundary(ptf, p, iF, mapper)
{}


NeumannCoupledBoundary::
NeumannCoupledBoundary
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
   CoupledBoundary(p, iF, dict)
{

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


NeumannCoupledBoundary::
NeumannCoupledBoundary
(
    const NeumannCoupledBoundary& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    CoupledBoundary(wtcsf, iF)
{}


NeumannCoupledBoundary::
NeumannCoupledBoundary
(
    const 
    NeumannCoupledBoundary& wtcsf
)
:
    CoupledBoundary(wtcsf)
{}


void NeumannCoupledBoundary::updateCoeffs()
{
    if (updated() || this->size() == 0)
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType(); 
   
    // Temperature field at boundary
    const scalarField& Tp = *this;
    
    // get cell centre coordinates of boundary cells
    vectorField coords = patch().Cf();


    // push temperature field at boundary to neighbour
    forAll(Tp, faceI){
        push(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), Tp[faceI], "temp");
    }  

    commit();


    // Thermal conductivity field at boundary
    const scalarField kappaTp = kappa(Tp);


    // fetch heat flux at boundary from neighbour
    scalarField nbrFluxFld = scalarField(this->size());
    forAll(nbrFluxFld, faceI){
        nbrFluxFld[faceI] = fetch(coords[faceI].x(), coords[faceI].y(), coords[faceI].z(), "flux");
    }


    // set reference gradient
    this->refGrad() = -nbrFluxFld/kappaTp;
    this->valueFraction() = 0;

    CoupledBoundary::updateCoeffs();
    UPstream::msgType(oldTag);  // Restore tag
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    
NeumannCoupledBoundary
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


// ************************************************************************* //
