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

  \*---------------------------------------------------------------------------*/

#include "timeVaryingHeatFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "Time.H"
#include "AverageIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  // declare specialization within 'Foam' namespace
  template<>
  const char* NamedEnum
  <
    Foam::incompressible::
    timeVaryingHeatFluxFvPatchScalarField::heatSourceType,
    2
    >::names[] =
    {
      "power",
      "flux"
    };
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

  namespace incompressible
  {

    // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    const NamedEnum
    <
      timeVaryingHeatFluxFvPatchScalarField::heatSourceType,
      2
      > timeVaryingHeatFluxFvPatchScalarField::heatSourceTypeNames_;


    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    timeVaryingHeatFluxFvPatchScalarField::
    timeVaryingHeatFluxFvPatchScalarField
    (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF
     )
      :
      fixedGradientFvPatchScalarField(p, iF),
      heatSource_(hsPower),
      q_(p.size(), 0.0),
      perturb_(0),
      Prt(0.7),
      fieldTableName_(iF.name()),
      radfluxName_("T"), 
      mapperPtr_(NULL),
      sampleTimes_(0),
      startSampleTime_(-1),
      startSampledValues_(0),
      endSampleTime_(-1),
      endSampledValues_(0),
      startFluxValues_(0),
      endFluxValues_(0) 
    {}


    timeVaryingHeatFluxFvPatchScalarField::
    timeVaryingHeatFluxFvPatchScalarField
    (
     const timeVaryingHeatFluxFvPatchScalarField& ptf,
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const fvPatchFieldMapper& mapper
     )
      :
      fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
      heatSource_(ptf.heatSource_),
      q_(ptf.q_, mapper),
      perturb_(ptf.perturb_),
      Prt(ptf.Prt),
      fieldTableName_(ptf.fieldTableName_),
      radfluxName_("T"), 
      mapMethod_(ptf.mapMethod_),
      mapperPtr_(NULL),
      sampleTimes_(0),
      startSampleTime_(-1),
      startSampledValues_(0),
      endSampleTime_(-1),
      endSampledValues_(0),
      startFluxValues_(0),
      endFluxValues_(0)

    {}


    timeVaryingHeatFluxFvPatchScalarField::
    timeVaryingHeatFluxFvPatchScalarField
    (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict
     )
      :
      fixedGradientFvPatchScalarField(p, iF),
      heatSource_(heatSourceTypeNames_.read(dict.lookup("heatSource"))),
      q_("q", dict, p.size()),
      perturb_(dict.lookupOrDefault("perturb", 1e-5)),
      Prt(dict.lookupOrDefault("Prt", 0.7)),
      fieldTableName_(iF.name()),
      radfluxName_("T"), 
      mapMethod_
      (
       dict.lookupOrDefault<word>
       (
	"mapMethod",
	"planarInterpolation"
        )
       ),
      mapperPtr_(NULL),
      sampleTimes_(0),
      startSampleTime_(-1),
      startSampledValues_(0),
      endSampleTime_(-1),
      endSampledValues_(0),
      startFluxValues_(0),
      endFluxValues_(0)

    {
      if
	(
	 mapMethod_ != "planarInterpolation"
	 && mapMethod_ != "nearest"
	 )
	{
	  FatalIOErrorIn
	    (
	     "timeVaryingMappedFixedValueFvPatchField<Type>::\n"
	     "timeVaryingMappedFixedValueFvPatchField\n"
	     "(\n"
	     "    const fvPatch&\n"
	     "    const DimensionedField<Type, volMesh>&\n"
	     "    const dictionary&\n"
	     ")\n",
	     dict
	     )   << "mapMethod should be one of 'planarInterpolation'"
		 << ", 'nearest'" << exit(FatalIOError);
	}
      if (dict.found("value") && dict.found("gradient"))
	{
	  fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
	  gradient() = Field<scalar>("gradient", dict, p.size());
	}
      else
	{
	  fvPatchField<scalar>::operator=(patchInternalField());
	  gradient() = 0.0;
	}
    

    }


    timeVaryingHeatFluxFvPatchScalarField::
    timeVaryingHeatFluxFvPatchScalarField
    (
     const timeVaryingHeatFluxFvPatchScalarField& thftpsf
     )
      :
      fixedGradientFvPatchScalarField(thftpsf),
      heatSource_(thftpsf.heatSource_),
      q_(thftpsf.q_),
      perturb_(thftpsf.perturb_),    
      Prt(thftpsf.Prt),
      fieldTableName_(thftpsf.fieldTableName_),
      radfluxName_("T"), 
      mapMethod_(thftpsf.mapMethod_),
      mapperPtr_(NULL),
      sampleTimes_(thftpsf.sampleTimes_),
      startSampleTime_(thftpsf.startSampleTime_),
      startSampledValues_(thftpsf.startSampledValues_),
      endSampleTime_(thftpsf.endSampleTime_),
      endSampledValues_(thftpsf.endSampledValues_),
      startFluxValues_(thftpsf.startFluxValues_),
      endFluxValues_(thftpsf.endFluxValues_)
  
    {}


    timeVaryingHeatFluxFvPatchScalarField::
    timeVaryingHeatFluxFvPatchScalarField
    (
     const timeVaryingHeatFluxFvPatchScalarField& thftpsf,
     const DimensionedField<scalar, volMesh>& iF
     )
      :
      fixedGradientFvPatchScalarField(thftpsf, iF),
      heatSource_(thftpsf.heatSource_),
      q_(thftpsf.q_),
      perturb_(thftpsf.perturb_),    
      Prt(thftpsf.Prt),
      fieldTableName_(thftpsf.fieldTableName_),
      radfluxName_("T"), 
      mapMethod_(thftpsf.mapMethod_),
      mapperPtr_(NULL),
      sampleTimes_(thftpsf.sampleTimes_),
      startSampleTime_(thftpsf.startSampleTime_),
      startSampledValues_(thftpsf.startSampledValues_),
      endSampleTime_(thftpsf.endSampleTime_),
      endSampledValues_(thftpsf.endSampledValues_),
      startFluxValues_(thftpsf.startFluxValues_),
      endFluxValues_(thftpsf.endFluxValues_)
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void timeVaryingHeatFluxFvPatchScalarField::autoMap
    (
     const fvPatchFieldMapper& m
     )
    {
      scalarField::autoMap(m);
      q_.autoMap(m);   
      if (startSampledValues_.size())
	{
	  startSampledValues_.autoMap(m);
	  endSampledValues_.autoMap(m);
          startFluxValues_.autoMap(m);
          endFluxValues_.autoMap(m);
	}
      // Clear interpolator
      mapperPtr_.clear();
      startSampleTime_ = -1;
      endSampleTime_ = -1;
    }


    void timeVaryingHeatFluxFvPatchScalarField::rmap
    (
     const fvPatchScalarField& ptf,
     const labelList& addr
     )
    {
      fixedGradientFvPatchScalarField::rmap(ptf, addr);

      const timeVaryingHeatFluxFvPatchScalarField& thftptf =
        refCast<const timeVaryingHeatFluxFvPatchScalarField>
        (
	 ptf
	 );

      q_.rmap(thftptf.q_, addr);
      startSampledValues_.rmap(thftptf.startSampledValues_, addr);
      endSampledValues_.rmap(thftptf.endSampledValues_, addr);
      startFluxValues_.rmap(thftptf.startFluxValues_, addr);
      endFluxValues_.rmap(thftptf.endFluxValues_, addr);

      // Clear interpolator
      mapperPtr_.clear();
      startSampleTime_ = -1;
      endSampleTime_ = -1;    
    }

    void timeVaryingHeatFluxFvPatchScalarField::checkTable()
    {
      // Initialise
      if (mapperPtr_.empty())
	{
	  pointIOField samplePoints
	    (
	     IOobject
	     (
	      "points",
	      this->db().time().constant(),
	      "boundaryData"/this->patch().name(),
	      this->db(),
	      IOobject::MUST_READ,
	      IOobject::AUTO_WRITE,
	      false
	      )
	     );

	  const fileName samplePointsFile = samplePoints.filePath();
	  if (debug)
	    {
	      Info<< "timeVaryingMappedFixedValueFvPatchField :"
		  << " Read " << samplePoints.size() << " sample points from "
		  << samplePointsFile << endl;
	    }
	  if (perturb_ < 1e-5)
	    {
	      Info<<"Overwriting the perturbation as it is small"<<endl;
	      perturb_=1e-5;
	    }
	  // Allocate the interpolator
	  mapperPtr_.reset
	    (
	     new pointToPointPlanarInterpolation
	     (
	      samplePoints,
	      this->patch().patch().faceCentres(),
	      perturb_
	      )
	     );

	  // Read the times for which data is available
	  const fileName samplePointsDir = samplePointsFile.path();
	  sampleTimes_ = Time::findTimes(samplePointsDir);
	  if (debug)
	    {
	      Info<< "timeVaryingMappedFixedValueFvPatchField : In directory "
		  << samplePointsDir << " found times "
		  << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
		  << endl;
	    }        
	}


      // Find current time in sampleTimes
      label lo = -1;
      label hi = -1;

      bool foundTime = mapperPtr_().findTime
	(
	 sampleTimes_,
	 startSampleTime_,
	 this->db().time().value(),
	 lo,
	 hi
	 );

      if (!foundTime)
	{
	  FatalErrorIn
	    (
	     "timeVaryingMappedFixedValueFvPatchField<Type>::checkTable()"
	     )   << "Cannot find starting sampling values for current time "
		 << this->db().time().value() << nl
		 << "Have sampling values for times "
		 << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
		 << "In directory "
		 <<  this->db().time().constant()/"boundaryData"/this->patch().name()
		 << "\n    on patch " << this->patch().name()
		 << " of field " << fieldTableName_
		 << exit(FatalError);
	}


      // Update sampled data fields.

      if (lo != startSampleTime_)
	{
	  startSampleTime_ = lo;

	  if (startSampleTime_ == endSampleTime_)
	    {
	      if (debug)
		{
		  Pout<< "checkTable : Setting startValues to (already read) "
		      <<   "boundaryData"
		    /this->patch().name()
		    /sampleTimes_[startSampleTime_].name()
		      << endl;
		}
	      // No need to reread since are end values
	      startSampledValues_ = endSampledValues_;
	    }
	  else
	    {
	      if (debug)
		{
		  Pout<< "checkTable : Reading startValues from "
		      <<   "boundaryData"
		    /this->patch().name()
		    /sampleTimes_[lo].name()
		      << endl;
		}
	      // Reread values and interpolate
	      AverageIOField<scalar> vals
		(
		 IOobject
		 (
		  fieldTableName_,
		  this->db().time().constant(),
		  "boundaryData"
		  /this->patch().name()
		  /sampleTimes_[startSampleTime_].name(),
		  this->db(),
		  IOobject::MUST_READ,
		  IOobject::AUTO_WRITE,
		  false
		  )
		 );

	      AverageIOField<scalar> fluxvals
		(
                 IOobject
                 (
                  radfluxName_,
                  this->db().time().constant(),
                  "boundaryData"
                  /this->patch().name()
                  /sampleTimes_[startSampleTime_].name(),
                  this->db(),
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE,
                  false
                  )
		 );

	      if (vals.size() != mapperPtr_().sourceSize())
		{
		  FatalErrorIn
		    (
		     "timeVaryingMappedFixedValueFvPatchField<Type>::"
		     "checkTable()"
		     )   << "Number of values (" << vals.size()
			 << ") differs from the number of points ("
			 <<  mapperPtr_().sourceSize()
			 << ") in file " << vals.objectPath() << exit(FatalError);
		}
              if (fluxvals.size() != mapperPtr_().sourceSize())
                {
                  FatalErrorIn
                    (
                     "timeVaryingMappedFixedValueFvPatchField<Type>::"
                     "checkTable()"
                     )   << "Number of values (" << vals.size()
                         << ") differs from the number of points ("
                         <<  mapperPtr_().sourceSize()
                         << ") in file " << vals.objectPath() << exit(FatalError);
                }

	      startSampledValues_ = mapperPtr_().interpolate(vals);
	      startFluxValues_ = mapperPtr_().interpolate(fluxvals);
	    }
	}

      if (hi != endSampleTime_)
	{
	  endSampleTime_ = hi;

	  if (endSampleTime_ == -1)
	    {
	      if (debug)
		{
		  Pout<< "checkTable : Clearing endValues" << endl;
		}
	      endSampledValues_.clear();
	      endFluxValues_.clear();
	    }
	  else
	    {
	      if (debug)
		{
		  Pout<< "checkTable : Reading endValues from "
		      <<   "boundaryData"
		    /this->patch().name()
		    /sampleTimes_[endSampleTime_].name()
		      << endl;
		}
	      // Reread values and interpolate
	      AverageIOField<scalar> vals
		(
		 IOobject
		 (
		  fieldTableName_,
		  this->db().time().constant(),
		  "boundaryData"
		  /this->patch().name()
		  /sampleTimes_[endSampleTime_].name(),
		  this->db(),
		  IOobject::MUST_READ,
		  IOobject::AUTO_WRITE,
		  false
		  )
		 );

              AverageIOField<scalar> fluxvals
                (
                 IOobject
                 (
                  radfluxName_,
                  this->db().time().constant(),
                  "boundaryData"
                  /this->patch().name()
                  /sampleTimes_[startSampleTime_].name(),
                  this->db(),
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE,
                  false
                  )
                 );

	      if (vals.size() != mapperPtr_().sourceSize())
		{
		  FatalErrorIn
		    (
		     "timeVaryingMappedFixedValueFvPatchField<Type>::"
		     "checkTable()"
		     )   << "Number of values (" << vals.size()
			 << ") differs from the number of points ("
			 <<  mapperPtr_().sourceSize()
			 << ") in file " << vals.objectPath() << exit(FatalError);
		}
              if (fluxvals.size() != mapperPtr_().sourceSize())
                {
                  FatalErrorIn
                    (
                     "timeVaryingMappedFixedValueFvPatchField<Type>::"
                     "checkTable()"
                     )   << "Number of values (" << vals.size()
                         << ") differs from the number of points ("
                         <<  mapperPtr_().sourceSize()
                         << ") in file " << vals.objectPath() << exit(FatalError);
                }

	      endSampledValues_ = mapperPtr_().interpolate(vals);
	      endFluxValues_ = mapperPtr_().interpolate(fluxvals);
	    }
	}
    }



    void timeVaryingHeatFluxFvPatchScalarField::updateCoeffs()
    {
      if (updated())
	{
	  return;
	}
      checkTable();
      const label patchI = this->patch().index();
      const turbulenceModel& turbulence =db().lookupObject<turbulenceModel>("turbulenceModel");
      const tmp<volScalarField> tnut = turbulence.nut();
      const volScalarField& nut = tnut();
      const scalarField& nutw = nut.boundaryField()[patchI];
      //const volScalarField blanking_=db().lookupObject<volScalarField>("treeblanking"); 
      //const scalarField blankingw_=blanking_.boundaryField()[patchI].patchInternalField();       
      //const volScalarField sb_=db().lookupObject<volScalarField>("solarblanking"); 
      //const scalarField sbw_=sb_.boundaryField()[patchI].patchInternalField();       
      // retrieve (constant) specific heat capacity from transport dictionary
      const IOdictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
      const scalar rhoCp0(readScalar(transportProperties.lookup(this->patch().name() & ".rhoCp0")));    
      if (endSampleTime_ == -1)
	{
	  if (debug)
	    {
	      Pout<< "updateCoeffs : Sampled, non-interpolated values"
		  << " from start time:"
		  << sampleTimes_[startSampleTime_].name() << nl;
	    }
	  //	  q_=(startSampledValues_);
	  q_=(startFluxValues_);

	}
      else
	{
	  scalar start = sampleTimes_[startSampleTime_].value();
	  scalar end = sampleTimes_[endSampleTime_].value();
	  scalar s = (this->db().time().value() - start)/(end - start);
	  //	  q_=(1 - s)*startSampledValues_ + s*endSampledValues_;  
	  q_=(1 - s)*startFluxValues_ + s*endFluxValues_;  
	  if (debug)
	    {
	      Pout<< "updateCoeffs : Sampled, interpolated values"
		  << " between start time:"
		  << sampleTimes_[startSampleTime_].name()
		  << " and end time:" << sampleTimes_[endSampleTime_].name()
		  << " with weight:" << s << endl;
	    }
        
	}    

      switch (heatSource_)
	{
		      
        case hsPower:
	  {
            const scalar Ap = gSum(patch().magSf());
            gradient() = Prt*q_/(Ap*rhoCp0*(nutw+5e-5));
            break;
	  }
        case hsFlux:
	  {
	    gradient() = Prt*q_/(rhoCp0*(nutw+5e-5));
	    //	  Pout<<"Minimum nut:"<<this->patch().name()<<"  "<<min(nutw)<<endl;
	    break;
	  }
        default:
	  {
            FatalErrorIn
	      (
	       "timeVaryingHeatFluxFvPatchScalarField"
	       "("
	       "const fvPatch&, "
	       "const DimensionedField<scalar, volMesh>&, "
	       "const dictionary&"
	       ")"
	       )   << "Unknown heat source type. Valid types are: "
		   << heatSourceTypeNames_ << nl << exit(FatalError);
	  }
	}
      if (debug)
	{
	  Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
	      << " max:" << gMax(*this) << endl;
	}
      fixedGradientFvPatchScalarField::updateCoeffs();
    }


    void timeVaryingHeatFluxFvPatchScalarField::write(Ostream& os) const
    {
      fixedGradientFvPatchScalarField::write(os);
      os.writeKeyword("heatSource") << heatSourceTypeNames_[heatSource_]
				    << token::END_STATEMENT << nl;
      os.writeKeyword("Prt") << Prt << token::END_STATEMENT << nl;
      q_.writeEntry("q", os);
      writeEntry("value", os);
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
     fvPatchScalarField,
     timeVaryingHeatFluxFvPatchScalarField
     );


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  } // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
