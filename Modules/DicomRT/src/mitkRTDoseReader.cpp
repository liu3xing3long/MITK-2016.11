/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


#include "mitkRTDoseReader.h"

#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkDICOMFileReaderSelector.h>
#include <mitkDICOMFileReader.h>
#include <mitkRTConstants.h>
#include <mitkImageStatisticsHolder.h>
#include <mitkIOMimeTypes.h>
#include <mitkDICOMTagPath.h>

#include "usModuleContext.h"
#include "usGetModuleContext.h"

#include "dcmtk/dcmrt/drtdose.h"

#include <itkShiftScaleImageFilter.h>
#include <itkCastImageFilter.h>

namespace mitk
{

    RTDoseReader::RTDoseReader() : AbstractFileReader(IOMimeTypes::DICOM_MIMETYPE_NAME(), "DICOM RTDose File Reader") {
        m_ServiceReg = RegisterService();
    }

    RTDoseReader::RTDoseReader(const RTDoseReader& other) : mitk::AbstractFileReader(other)
    {

    }

    RTDoseReader::~RTDoseReader(){}

    template<typename TPixel, unsigned int VImageDimension>
    void RTDoseReader::MultiplyGridScaling(itk::Image<TPixel, VImageDimension>* image, float gridscale)
    {
        typedef itk::Image<Float32, VImageDimension> OutputImageType;
        typedef itk::Image<TPixel, VImageDimension> InputImageType;

        typedef itk::CastImageFilter<InputImageType, OutputImageType> CastFilterType;
        typedef itk::ShiftScaleImageFilter<OutputImageType, OutputImageType> ScaleFilterType;
        typename CastFilterType::Pointer castFilter = CastFilterType::New();
        typename ScaleFilterType::Pointer scaleFilter = ScaleFilterType::New();

        castFilter->SetInput(image);
        scaleFilter->SetInput(castFilter->GetOutput());
        scaleFilter->SetScale(gridscale);
        scaleFilter->Update();
        typename OutputImageType::Pointer scaledOutput = scaleFilter->GetOutput();
        this->scaledDoseImage = mitk::Image::New();

        mitk::CastToMitkImage(scaledOutput, this->scaledDoseImage);
    }

    mitk::IDICOMTagsOfInterest* RTDoseReader::GetDicomTagsOfInterestService()
    {
        mitk::IDICOMTagsOfInterest* result = nullptr;

        std::vector<us::ServiceReference<mitk::IDICOMTagsOfInterest> > toiRegisters = us::GetModuleContext()->GetServiceReferences<mitk::IDICOMTagsOfInterest>();
        if (!toiRegisters.empty())
        {
            if (toiRegisters.size() > 1)
            {
                MITK_WARN << "Multiple DICOM tags of interest services found. Using just one.";
            }
            result = us::GetModuleContext()->GetService<mitk::IDICOMTagsOfInterest>(toiRegisters.front());
        }

        return result;
    }

    std::vector<itk::SmartPointer<BaseData> > RTDoseReader::Read()
    {
        std::vector<itk::SmartPointer<mitk::BaseData> > result;

        DICOMTag referencedRTPlan(0x300c, 0x0002);
        mitk::IDICOMTagsOfInterest* toiSrv = GetDicomTagsOfInterestService();
        if (toiSrv)
        {
            toiSrv->AddTagOfInterest(referencedRTPlan);
        }

        std::string location = GetInputLocation();
        mitk::DICOMFileReaderSelector::Pointer selector = mitk::DICOMFileReaderSelector::New();
        selector->LoadBuiltIn3DConfigs();
        selector->SetInputFiles({ location });

        mitk::DICOMFileReader::Pointer reader = selector->GetFirstReaderWithMinimumNumberOfOutputImages();
        reader->SetAdditionalTagsOfInterest(toiSrv->GetTagsOfInterest());

        reader->SetInputFiles({ location });
        reader->AnalyzeInputFiles();
        reader->LoadImages();

        if (reader->GetNumberOfOutputs() == 0){
            MITK_ERROR << "Could not determine a DICOM reader for this file" << std::endl;
            return result;
        }

        const mitk::DICOMImageBlockDescriptor& desc = reader->GetOutput(0);

        mitk::Image::Pointer originalImage = desc.GetMitkImage();

        if (originalImage.IsNull())
        {
            MITK_ERROR << "Error reading the RTDOSE file in mitk::DicomFileReader" << std::endl;
            return result;
        }

        DcmFileFormat fileformat;
        OFCondition outp = fileformat.loadFile(location.c_str(), EXS_Unknown);
        if (outp.bad())
        {
            MITK_ERROR << "Error reading the RTDOSE file in DCMTK" << std::endl;
            return result;
        }
        DcmDataset *dataset = fileformat.getDataset();

        DRTDoseIOD doseObject;
        OFCondition DCMTKresult = doseObject.read(*dataset);

        if (DCMTKresult.bad())
        {
            MITK_ERROR << "Error reading the RTDOSE file in DCMTK" << std::endl;
            return result;
        }

        OFString gridScaling;
        Float32 gridscale;

        doseObject.getDoseGridScaling(gridScaling);
        gridscale = OFStandard::atof(gridScaling.c_str());

        AccessByItk_1(originalImage, MultiplyGridScaling, gridscale);

        auto statistics = this->scaledDoseImage->GetStatistics();
        double maxDose = statistics->GetScalarValueMax();

        this->scaledDoseImage->SetPropertyList(originalImage->GetPropertyList());
        this->scaledDoseImage->SetProperty(mitk::RTConstants::PRESCRIBED_DOSE_PROPERTY_NAME.c_str(), mitk::GenericProperty<double>::New(0.8*maxDose));

        result.push_back(this->scaledDoseImage.GetPointer());
        return result;
    }

    RTDoseReader* RTDoseReader::Clone() const
    {
        return new RTDoseReader(*this);
    }

}
