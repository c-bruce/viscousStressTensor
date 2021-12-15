/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::viscousStressTensor

Group
    grpFieldFunctionObjects

Description
    Computes the viscousStressTensor for the porousIcoFoam solver.

    Returns -rhoRef_*nu*dev(twoSymm(fvc::grad(U))) tensor field.

Usage
    Minimal example of system/viscousStressTensor.:
    \verbatim
    type		viscousStressTensor;
    libs		("fieldFunctionObjects" "libmyViscousStrainRateFunctionObject.so");

    field		U;
    rho		    1;

    executeControl	writeTime;
    writeControl	writeTime;
    \endverbatim

    Note: the "viscousStressTensor" file must be placed in the system directory to be
    used with the postProcess utility.

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func viscousStressTensor
    \endverbatim    

SourceFiles
    viscousStressTensor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_viscousStressTensor_H
#define functionObjects_viscousStressTensor_H

#include "fieldExpression.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class viscousStressTensor Declaration
\*---------------------------------------------------------------------------*/

class viscousStressTensor
:
    public fieldExpression,
    public logFiles
{

private:
    // Read from dictionary

        //- Reference density for incompressible calculations
        scalar rhoRef_;

    // Private Member Functions

        using logFiles::files; // OF2106

        //- Calculate the viscousStressTensor field and return true if successful
        virtual bool calc();

        //- Return the effective viscous stress
        tmp<volSymmTensorField> devTau() const;

public:

    //- Runtime type information
    TypeName("viscousStressTensor");

    // Constructors

        //- Construct from Time and dictionary
        viscousStressTensor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        viscousStressTensor(const viscousStressTensor&) = delete;

        //- No copy assignment
        void operator=(const viscousStressTensor&) = delete;

    //- Destructor
    virtual ~viscousStressTensor() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //