/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "viscousStressTensor.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(viscousStressTensor, 0);
    addToRunTimeSelectionTable(functionObject, viscousStressTensor, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::viscousStressTensor::devTau() const
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "constant/transportProperties",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        transportProperties // OF2106
    );

    const volVectorField& U = obr_.lookupObject<volVectorField>(fieldName_);

    return -rhoRef_*nu*dev(twoSymm(fvc::grad(U)));
}

bool Foam::functionObjects::viscousStressTensor::calc()
{
    // Get viscous stress tensor field
    tmp<volSymmTensorField> tdevTau = devTau();

    return store
    (
        resultName_,
        tdevTau
    );

    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::viscousStressTensor::viscousStressTensor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U"),
    logFiles(obr_, name),
    rhoRef_(dict.get<scalar>("rho"))
{
    setResultName(typeName, fieldName_);
}

// ************************************************************************* //