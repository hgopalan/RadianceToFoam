/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
  \\/     M anipulation  |
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
*/

/*
  !---------------------------------------------------------------------------!  
  ! Time Varying Version of Surface Temperature Flux Boundary Condition       !
  !                                                                           !
  ! Description                                                               !
  ! This boundary condition merges the time varying mapped fixedValue BC and  !
  ! surface temperature heat-flux BC to provide a time-varying heat flux BC.  !
  ! Rest of the details are very similar to the surfaceTurbulentHeatFlux BC.  !
  !                                                                           !
  !                                                                           !
  ! History                                                                   !
  ! 1. October 2015: Created by merging timeVaryingmappedFixedValue and       !
  !                  surface heat flux boundary condition.                    !
  ! 2. May 2016: There was a bug in array indexing. Fixed now.                !
  !                                                                           !
  !                                                                           !
  ! Notes                                                                     !
  ! 1. The constant/boundaryData/0.. and other time-instant folders should    !
  !    contain the heat-flux as the input values with filename T. This issue  !
  !    will be updated if required.                                           !    
  ! 2. Use the accompanying reader.py utility to create the files from the    !
  !    output provided from Radiance software. Currently, the script writes   !
  !    out information for only one time-step. Modify, if necessary.          !
  ! 3. Only the implementation has been verified. Please test it thoroughly   !
  !    and let me know of any issues.                                         !
  ! 4. The density*specific heat capacity rhoCp0 is read from the transport   !
  !    dictionary entry. Be consistent with units.                            ! 
  !                                                                           !   
  !                                                                           !
  ! Authors                                                                   !
  ! Harish Gopalan (gopalanh@ihpc.a-star.edu.sg)                              !
  !---------------------------------------------------------------------------!
*/

/* Specifying the Boundary Condition 
   buildingWall
   {
   type            timeVaryingHeatFlux;
   heatSource flux; //flux or power
   q uniform 20.0; // Dummy place-holder for non-uniform heat flux. Over-written automatically. 
   Prt 0.7; // Turbulent Prandtl number  
   // mapMethod specifies interpolation. Interpolation uses delanuay triangulation. 
   mapMethod nearest; // nearest or planarInterpolation 
   value uniform 2; // Dummy place-holder for actal temperature gradient. Over-written automatically. 
   }
*/



#ifndef timeVaryingHeatFluxFvPatchScalarFields_H
#define timeVaryingHeatFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "NamedEnum.H"
#include "FixedList.H"
#include "instantList.H"
#include "pointToPointPlanarInterpolation.H"
#include "DataEntry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  namespace incompressible
  {

    /*---------------------------------------------------------------------------*\
      Class timeVaryingHeatFluxFvPatchScalarField Declaration
      \*---------------------------------------------------------------------------*/

    class timeVaryingHeatFluxFvPatchScalarField
      :
      public fixedGradientFvPatchScalarField
    {
    public:

      // Data types

      //- Enumeration listing the possible hest source input modes
      enum heatSourceType
        {
	  hsPower,
	  hsFlux
        };


    private:

      // Private data

      //- Heat source type names
      static const NamedEnum<heatSourceType, 2> heatSourceTypeNames_;

      //- Heat source type
      heatSourceType heatSource_;

      //- Heat power [W] or flux [W/m2]
      scalarField q_;

      //[HG] Perturb the mesh to avoid breakdown of Delanuay triangulation 
      scalar perturb_;

      // Turbulent Prandtl number 
      scalar Prt;

      //- Name of the field data table, defaults to the name of the field.
      // [HG] *************** Need to modify it to read the heat-flux filename in future. **************
      word fieldTableName_;

      // HG - Added on May 13th 2016
      word radfluxName_;
      // End - Addition 

      //- Interpolation scheme to use
      word mapMethod_;
      
      //- 2D interpolation pointer
      autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

      //- List of boundaryData time directories
      instantList sampleTimes_;

      //- Current starting index in sampleTimes
      label startSampleTime_;

      //- Interpolated values from startSampleTime
      scalarField startSampledValues_;

      //- Current end index in sampleTimes
      label endSampleTime_;

      //- Interpolated values from endSampleTime
      scalarField endSampledValues_;

      // HG - Added on May 13th 2016
      scalarField startFluxValues_;
      scalarField endFluxValues_;

    public:

      //- Runtime type information
      TypeName("timeVaryingHeatFlux");


      // Constructors

      //- Construct from patch and internal field
      timeVaryingHeatFluxFvPatchScalarField
      (
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&
       );

      //- Construct from patch, internal field and dictionary
      timeVaryingHeatFluxFvPatchScalarField
      (
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&,
       const dictionary&
       );

      //- Construct by mapping given
      //  timeVaryingHeatFluxFvPatchScalarField onto
      //  a new patch
      timeVaryingHeatFluxFvPatchScalarField
      (
       const timeVaryingHeatFluxFvPatchScalarField&,
       const fvPatch&,
       const DimensionedField<scalar, volMesh>&,
       const fvPatchFieldMapper&
       );

      //- Construct as copy
      timeVaryingHeatFluxFvPatchScalarField
      (
       const timeVaryingHeatFluxFvPatchScalarField&
       );

      //- Construct and return a clone
      virtual tmp<fvPatchScalarField> clone() const
      {
	return tmp<fvPatchScalarField>
	  (
	   new timeVaryingHeatFluxFvPatchScalarField(*this)
	   );
      }

      //- Construct as copy setting internal field reference
      timeVaryingHeatFluxFvPatchScalarField
      (
       const timeVaryingHeatFluxFvPatchScalarField&,
       const DimensionedField<scalar, volMesh>&
       );

      //- Construct and return a clone setting internal field reference
      virtual tmp<fvPatchScalarField> clone
      (
       const DimensionedField<scalar, volMesh>& iF
       ) const
      {
	return tmp<fvPatchScalarField>
	  (
	   new timeVaryingHeatFluxFvPatchScalarField
	   (
	    *this,
	    iF
	    )
	   );
      }


      // Member functions

      // Mapping functions

      //- Map (and resize as needed) from self given a mapping object
      virtual void autoMap(const fvPatchFieldMapper&);

      //- Reverse map the given fvPatchField onto this fvPatchField
      virtual void rmap
      (
       const fvPatchScalarField&,
       const labelList&
       );


      // Evaluation functions
      void checkTable();

      //- Update the coefficients associated with the patch field
      virtual void updateCoeffs();


      // I-O

      //- Write
      virtual void write(Ostream&) const;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
