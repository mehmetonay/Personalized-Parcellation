
//ITK LIBRARIES
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkVariableLengthVector.h>

//STANDARD LIBRARIES
#include <iostream>
#include <fstream>
#include <string>

#include "Parc.h"

int main(int argc, char *argv[])
{
//    Type definitions for parcel and fmri images
    typedef itk::Image<unsigned int,3> ParcelImageType;
    typedef itk::Image<float,4> FMRIImageType;
    
//    Type definitions of readers for parcel and fmri images
    typedef itk::ImageFileReader<ParcelImageType> ParcelImageReaderType;
    typedef itk::ImageFileReader<FMRIImageType> FMRIImageReaderType;
    
    std::vector<std::string> filenames_AD(13);
	filenames_AD[0] = "CAPA-00015";
	filenames_AD[1] = "CAPA-00037";
	filenames_AD[2] = "CAPA-00045";
	filenames_AD[3] = "CAPA-00050";
	filenames_AD[4] = "CAPA-00053";
	filenames_AD[5] = "CAPA-00058";
	filenames_AD[6] = "CAPA-00063";
	filenames_AD[7] = "CAPA-00073";
	filenames_AD[8] = "CAPA-00076";
	filenames_AD[9] = "CAPA-00077";
	filenames_AD[10] = "CAPA-00079";
	filenames_AD[11] = "CAPA-00091";
	filenames_AD[12] = "CAPA-00092";

	std::vector<std::string> filenames_Healthy(13);
	filenames_Healthy[0] = "CAPA-00004";
	filenames_Healthy[1] = "CAPA-00005";
	filenames_Healthy[2] = "CAPA-00018";
	filenames_Healthy[3] = "CAPA-00021";
	filenames_Healthy[4] = "CAPA-00022";
	filenames_Healthy[5] = "CAPA-00026";
	filenames_Healthy[6] = "CAPA-00034";
	filenames_Healthy[7] = "CAPA-00038";
	filenames_Healthy[8] = "CAPA-00040";
	filenames_Healthy[9] = "CAPA-00041";
	filenames_Healthy[10] = "CAPA-00044";
	filenames_Healthy[11] = "CAPA-00047";
	filenames_Healthy[12] = "CAPA-00049";
    
	for (int i = 0; i < filenames_AD.size(); i++)
	{
		std::cout << "LAMBDA  " << lamda << "  PATIENT: " << filenames_AD[i] << std::endl;
		//    Add paths for parcellation and fmri images
		//        std::string patientName = "CAPA-00004";
		std::string patientName = filenames_AD[i];

		//    Add paths for parcellation and fmri images
		std::string parcelImageFilePath = "/Users/Mehmet/Desktop/DATA/AD/" + patientName + "/Processed/" + patientName + ".par.nii";
		std::string fmriImageFilePath = "/Users/Mehmet/Desktop/DATA/AD/" + patientName + "/Processed/" + patientName + "_preprocessed_fmri.fmri.nii.gz";

		//    Parcel and fmri image pointers defined
		ParcelImageType::Pointer parcelImage = ParcelImageType::New();
		ParcelImageType::Pointer newParcelImage = ParcelImageType::New();
		FMRIImageType::Pointer fmriImage = FMRIImageType::New();

		//    Reading parcel image
		ParcelImageReaderType::Pointer parcelImageReader = ParcelImageReaderType::New();
		parcelImageReader->SetFileName(parcelImageFilePath);
		parcelImage = parcelImageReader->GetOutput();
		parcelImage->Update();
		std::cout << "\nParcel image with file path '" << parcelImageFilePath << "' is read.\n";
		std::cout << "Size of parcel image: " << parcelImage->GetLargestPossibleRegion().GetSize() << std::endl;

		//    Reading fmri image
		FMRIImageReaderType::Pointer fmriImageReader = FMRIImageReaderType::New();
		fmriImageReader->SetFileName(fmriImageFilePath);
		fmriImage = fmriImageReader->GetOutput();
		fmriImage->Update();
		std::cout << "\nFMRI image with file path '" << fmriImageFilePath << "' is read.\n";
		std::cout << "Size of fmri image: " << fmriImage->GetLargestPossibleRegion().GetSize() << std::endl;

		//    Creating a new Parcellation object
		Vav::Parcellation::Parcellation parcellation;
		parcellation.SetFMRIImage(fmriImage);
		parcellation.SetParcelImage(parcelImage);
		parcellation.SetBoundOfCortexLabels(11000, 13000);
		parcellation.Update();

		//yeniden eklemeler umut

		parcellation.SetMaxParcelPopulation(maxParcelPopulation); // In case the big parcels are to be divided
		parcellation.SetResultsFilename("/Users/Mehmet/Desktop/VAVLab/network/AD/svdlike/" + patientName + "_svdlike_L025_I30.txt");
		parcellation.doSomethingg();

		newParcelImage = parcellation.GetNewParcelImage();
		newParcelImage->Update();
		//      Writing new parcellation
		std::string outputFilename;
		outputFilename = "/Users/Mehmet/Desktop/VAVLab/network/AD/svdlike/" + patientName + "_svdlike_L025_I30.par.nii";
		typedef  itk::ImageFileWriter< ParcelImageType  > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(outputFilename);
		writer->SetInput(newParcelImage);
		writer->Update();
	}

    return EXIT_SUCCESS;
}
