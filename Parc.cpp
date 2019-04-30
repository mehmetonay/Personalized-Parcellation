#include "Parc.h"

Vav::Parcellation::Parcellation::Parcellation()
{

}

Vav::Parcellation::Parcellation::~Parcellation()
{

}



void Vav::Parcellation::Parcellation::SetParcelImage(ParcelImageType::Pointer p)
{
	/*
	After reading parcel image, it should be set to parcellation object using this function
	*/
	std::clock_t startTime = std::clock();
	rawParcelImage = p;
	rawParcelImage->Update();
	std::clock_t endTime = std::clock();
	std::cout << "Parcel image set in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";
}


void Vav::Parcellation::Parcellation::SetFMRIImage(RawFMRIImageType::Pointer p)
{
	/*
	After reading fmri image, it should be set to parcellation object using this function
	*/
	std::clock_t startTime = std::clock();
	rawFMRIImage = p;
	rawFMRIImage->Update();
	std::clock_t endTime = std::clock();
	std::cout << "FMRI image set in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";
}


void Vav::Parcellation::Parcellation::SetBoundOfCortexLabels(const unsigned int lower, const unsigned int upper)
{
	/*
	Parcellation image contains not only cortex of the human brain but also the other parts. 
	In order to specify which region the algorithm will work on, we need to specify the interval 
	that is used for cortex labelling. In our case, cortex is labelled with numbers between 11000 and 13000, 
	therefore we set the interval as 
		parcellation.SetBoundOfCortexLabels(11000, 13000);
	*/
	upperBoundOfCortexLabels = upper;
	lowerBoundOfCortexLabels = lower;
	std::cout << "Cortex label bounds set.\n";
}



bool IsLabelIncluded(Vav::Parcellation::Parcellation::ParcelImageType::PixelType a, std::vector<Vav::Parcellation::Parcellation::ParcelImageType::PixelType> &vec)
{
	/*
	This function is used to check if the labels vector, which contains all the cortex labels, contains a specific label.
	We used this function in preprocessing part when transforming cortex labels into numbers 1,2,3,... and set zero for the rest of the brain.
	*/
	if (std::find(vec.begin(), vec.end(), a) != vec.end()) {
		return true;
	}
	else {
		return false;
	}
}

Vav::Parcellation::Parcellation::ParcelImageType::PixelType GetNewLabel(Vav::Parcellation::Parcellation::ParcelImageType::PixelType oldLabel, std::vector<Vav::Parcellation::Parcellation::ParcelImageType::PixelType> &vec)
{
	/*
	We used this function while transforming labels of the parcel image into 1,2,3,etc for cortex and zero for the rest of the brain.
	*/
	if (oldLabel == 0) {
		return 0;
	}
	else {
		Vav::Parcellation::Parcellation::ParcelImageType::PixelType newLabel = 0, counter = 0;
		std::vector<Vav::Parcellation::Parcellation::ParcelImageType::PixelType>::iterator it;
		for (it = vec.begin(); it != vec.end(); it++) {
			++counter;
			if ((*it) == oldLabel) {
				newLabel = counter;
				break;
			}
		}
		return newLabel;
	}

}

void Vav::Parcellation::Parcellation::Update(){
	/*
	This function performs most of the preprocessing part. It 
	-	initializes images
	-	finds the labels used for cortex
	-	changes labels of the cortex to 1,2,3,...
	-	changes the type of the fMRI image
	*/

	//    Performance calculation variables declared
	std::clock_t startTime, endTime;

	std::cout << "Starting update...\n";

	//    Initialization of images
	startTime = std::clock();
	std::cout << "Starting initialization of images...\n";
	ParcelImageType::RegionType::IndexType parcelStartIndex;
	parcelStartIndex[0] = 0;
	parcelStartIndex[1] = 0;
	parcelStartIndex[2] = 0;
	ParcelImageType::RegionType::SizeType parcelSize;
	parcelSize = rawParcelImage->GetLargestPossibleRegion().GetSize();
	ParcelImageType::RegionType parcelRegion;
	parcelRegion.SetIndex(parcelStartIndex);
	parcelRegion.SetSize(parcelSize);
	initialParcelImage->SetRegions(parcelRegion);
	initialParcelImage->Allocate();
	parcelImage->SetRegions(parcelRegion);
	parcelImage->Allocate();
	fmriImage->SetRegions(parcelRegion);
	fmriImage->Allocate();
	DisplacementFieldType::RegionType::IndexType startIndexDisplacement;
	startIndexDisplacement.Fill(0);
	DisplacementFieldType::RegionType::SizeType sizeDisplacement;
	sizeDisplacement = parcelImage->GetLargestPossibleRegion().GetSize();
	DisplacementFieldType::RegionType regionDisplacement;
	regionDisplacement.SetSize(sizeDisplacement);
	regionDisplacement.SetIndex(startIndexDisplacement);
	displacementField->SetRegions(regionDisplacement);
	displacementField->Allocate();
	directionImage->SetRegions(regionDisplacement);
	directionImage->Allocate();
	endTime = std::clock();
	std::cout << "Images initialized in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";

	//    Labels of parcel image between lower and upper bounds is found
	startTime = std::clock();
	std::cout << "Starting searching labels between " << lowerBoundOfCortexLabels << " and " << upperBoundOfCortexLabels << " in the parcel image...\n";
	typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
	ParcelImageRegionIteratorType itRawParcel(rawParcelImage, rawParcelImage->GetLargestPossibleRegion());
	itRawParcel.GoToBegin();
	while (!itRawParcel.IsAtEnd()) {
		if (itRawParcel.Get() >= lowerBoundOfCortexLabels && itRawParcel.Get() <= upperBoundOfCortexLabels) {
			if (!IsLabelIncluded(itRawParcel.Get(), labels)) {
				labels.push_back(itRawParcel.Get());
			}
		}
		++itRawParcel;
	}
	numberOfSegments = labels.size();
	endTime = std::clock();
	std::cout << labels.size() << " labels found in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";

	//    Labels of parcel image is changed as 1,2,3,... and stored
	startTime = std::clock();
	std::cout << "Starting change of parcel image labels to 1,2,3,...\n";
	ParcelImageRegionIteratorType itInitialParcel(initialParcelImage, initialParcelImage->GetLargestPossibleRegion());
	itRawParcel.GoToBegin();
	itInitialParcel.GoToBegin();
	while (!itRawParcel.IsAtEnd()) {
		itInitialParcel.Set(GetNewLabel(itRawParcel.Get(), labels));
		++itRawParcel;
		++itInitialParcel;
	}

	ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
	itInitialParcel.GoToBegin();
	itParcel.GoToBegin();

	while (!itParcel.IsAtEnd()) 
	{
        itParcel.Set(itInitialParcel.Get());
        ++itParcel;
        ++itInitialParcel;
    }
	// parcelImage = initialParcelImage;



	endTime = std::clock();
	std::cout << "Labels of parcel images is changed as 1,2,3... and stored in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";


	//    4D FMRI image is transformed into a 3D vector image
	// We used vector multiplication to calculate correlation coefficients, for that, we need fmri
	// signals to be vectors. That is why we transform 4D fmri image into 3D image that contains time signals at each voxel.
	std::cout << "Starting 4D fmri to 3D fmri transformation...\n";
	startTime = clock();

	typedef itk::ImageLinearIteratorWithIndex<RawFMRIImageType> RawFMRILinearIteratorType;
	RawFMRILinearIteratorType itRawFMRI(rawFMRIImage, rawFMRIImage->GetLargestPossibleRegion());
	itRawFMRI.SetDirection(3);
	itRawFMRI.GoToBegin();

	typedef itk::ImageRegionIterator<FMRIImageType> FMRIRegionIteratorType;
	FMRIRegionIteratorType itFMRI(fmriImage, fmriImage->GetLargestPossibleRegion());
	itFMRI.GoToBegin();

	FMRIImageType::PixelType timeSignal;
	timeSignal.set_size(rawFMRIImage->GetLargestPossibleRegion().GetSize()[3]);
	FMRIImageType::PixelType::iterator itTimeSignal;

	while (!itFMRI.IsAtEnd()) {
		itRawFMRI.GoToBeginOfLine();
		itTimeSignal = timeSignal.begin();
		while (!itRawFMRI.IsAtEndOfLine()) {
			(*itTimeSignal) = itRawFMRI.Get();
			++itTimeSignal;
			++itRawFMRI;
		}
		itFMRI.Set(timeSignal);
		++itFMRI;
		itRawFMRI.NextLine();
	}

	//    Freeing memory for 4D fmri
	rawFMRIImage = NULL;
	endTime = clock();
	std::cout << "4D fmri to 3D fmri transformation completed in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";

	//    Starting subtracting mean from time signals...
	startTime = clock();
	itFMRI.GoToBegin();
	while (!itFMRI.IsAtEnd()) {
		itFMRI.Value().operator-=(itFMRI.Get().mean());
		//        Normalization of fmri signals below!!!!!!!!!!!!!!!!!!
		itFMRI.Value().normalize();
		++itFMRI;
	}
	endTime = clock();
	std::cout << "Means of fmri signals subtracted in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds.\n";
	std::cout << "Parcellation update done.\n";
}




std::vector<Vav::Parcellation::Parcellation::ParcelImageType::PixelType> Vav::Parcellation::Parcellation::GetLabels(){
	return labels;
}



Vav::Parcellation::Parcellation::ParcelImageType::Pointer Vav::Parcellation::Parcellation::GetInitialParcelImage(){
	return initialParcelImage;
}



Vav::Parcellation::Parcellation::FMRIImageType::Pointer Vav::Parcellation::Parcellation::GetVectorFMRIImage(){
	return fmriImage;
}


void Vav::Parcellation::Parcellation::createFinalParcelImage()
{
	 //Allocation of finalParcelImage

    ParcelImageType::RegionType::IndexType parcelStartIndex;
    parcelStartIndex[0] = 0;
    parcelStartIndex[1] = 0;
    parcelStartIndex[2] = 0;
    ParcelImageType::RegionType::SizeType parcelSize;
    parcelSize = rawParcelImage->GetLargestPossibleRegion().GetSize();
    ParcelImageType::RegionType parcelRegion;
    parcelRegion.SetIndex(parcelStartIndex);
    parcelRegion.SetSize(parcelSize);
    
    finalParcelImage->SetRegions(parcelRegion);
    finalParcelImage->SetSpacing(rawParcelImage->GetSpacing());
    finalParcelImage->SetOrigin(rawParcelImage->GetOrigin());
    finalParcelImage->SetDirection(rawParcelImage->GetDirection());
    finalParcelImage->Allocate();

//    Iterators defined
    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    ParcelImageRegionIteratorType itFinalParcel(finalParcelImage, finalParcelImage->GetLargestPossibleRegion());
    itFinalParcel.GoToBegin();
    ParcelImageRegionIteratorType itRawParcel(rawParcelImage, rawParcelImage->GetLargestPossibleRegion());
    itRawParcel.GoToBegin();

    while (!itRawParcel.IsAtEnd()) {
        itFinalParcel.Set(itRawParcel.Get());
        ++itRawParcel;
        ++itFinalParcel;
    }
    itRawParcel.GoToBegin();
    itFinalParcel.GoToBegin();

    while (!itParcel.IsAtEnd()) {
        if (itParcel.Get()>0) {
            itFinalParcel.Value() = labels[itParcel.Get() - 1];
        }
        ++itParcel;
        ++itFinalParcel;
    }
    
    // finalParcelImage->Update();
}


void Vav::Parcellation::Parcellation::findSegmentInformation()
{
	/*
	This function if used to generate segment representative signals, which is in this case, 
	summary signals obtained using SVD decomposition. In order to do this,
	*/
	std::cout << "Starting finding segment information..." << std::endl;
	std::clock_t startTime = std::clock();
	std::clock_t startTimeTop = std::clock();

	//    Initialization of numberOfSegmentVoxels and segmentMeans
	numberOfSegmentVoxels.set_size(numberOfSegments);
	numberOfSegmentVoxels.fill(0);
	FMRIImageType::PixelType initializationVector;
	FMRIImageType::IndexType index;
	index[0] = 0;
	index[1] = 0;
	index[2] = 0;
	initializationVector = fmriImage->GetPixel(index);
	initializationVector = fmriImage->GetPixel(index);
	initializationVector.operator*=(0);
	segmentMeans.clear();
	segmentSVDs.clear();

	for (unsigned int i = 0; i < numberOfSegments; i++)
	{
		segmentMeans.push_back(initializationVector);
		segmentSVDs.push_back(initializationVector);
	}

	std::clock_t endTime = std::clock();
	std::cout << "numberOfSegmentVoxels, segmentMeans and segmentSVDs containers initialized in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;

	//    Cumulating # of voxels and their fmri signals

	startTime = std::clock();
	
	typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
	ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
	itParcel.GoToBegin();
	typedef itk::ImageRegionIterator<FMRIImageType> FMRIImageRegionIterator;
	FMRIImageRegionIterator itFMRI(fmriImage, fmriImage->GetLargestPossibleRegion());
	itFMRI.GoToBegin();

	while (!itParcel.IsAtEnd())
	{
		if (itParcel.Get() > 0)
		{
			numberOfSegmentVoxels[itParcel.Get() - 1] += 1;
			segmentMeans[itParcel.Get() - 1].operator+=(itFMRI.Get());
		}
		++itParcel;
		++itFMRI;
	}

	endTime = std::clock();
	std::cout << "Related information cumulated in initialized containers in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;

	//    Finding mean signals using numberOfSegmentVoxels and segmentMeans cumulations
	startTime = std::clock();

	std::vector<FMRIImageType::PixelType>::iterator itSegmentMeans;
	itSegmentMeans = segmentMeans.begin();
	vnl_vector<int>::iterator itNumberOfSegmentVoxels;
	itNumberOfSegmentVoxels = numberOfSegmentVoxels.begin();

	while (itSegmentMeans != segmentMeans.end()) 
	{
		(*itSegmentMeans).operator/=((*itNumberOfSegmentVoxels));
		++itSegmentMeans;
		++itNumberOfSegmentVoxels;
	}

	endTime = std::clock();
	std::cout << "Dividing by # of voxels to obtain means performed in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;

	// Finding SVD summary signals of the parcels
	// All signals are put into a matrix as columns to find SVD decomposition
	startTime = std::clock();

	itParcel.GoToBegin();
	itFMRI.GoToBegin();

	typedef Eigen::MatrixXf MatrixType;
	typedef std::vector<vnl_vector<float>> _MatrixType;

	unsigned int f = 1;
	while (f < numberOfSegments+1)
	{
		_MatrixType x;
		x.clear();

		while (!itParcel.IsAtEnd())
		{
			if (itParcel.Get() == f)
			{
				x.push_back(itFMRI.Get());
			}
			++itParcel;
			++itFMRI;
		}
		
		// All signals are put into a matrix as columns to find SVD decomposition
		MatrixType _x(x[0].size(), x.size());
		
		for (int i = 0; i < x.size(); i++){
			for (int j = 0; j < x[0].size(); j++){
				_x(j, i) = x[i][j];
			}
		}

		

		Eigen::BDCSVD<MatrixType> SVD(_x, Eigen::DecompositionOptions::ComputeFullU);
		
		auto S = SVD.singularValues();
		int maxindex;
 		auto max = S.maxCoeff(&maxindex);
		auto vec = SVD.matrixU().col(maxindex);

		for (int i = 0; i < vec.size(); i++)
		{
			segmentSVDs[f - 1].put(i, vec(i));
		}

		// Due to the nature of SVD method here, left most singular vector might be found in the (-) direction.
		// However it needs to be positively correlated with the parcel mean signal. In such a case, basic sign change of the vector will work.
		if(dot_product(segmentSVDs[f - 1],segmentMeans[f - 1]) < 0)
		{
			segmentSVDs[f-1].operator*=(-1);
		}
		// summary signal found above is used to find the "intrinsic summary signal" which corresponds to the voxel signal with
		//the highest correlation with summary signal. It gives similar result in terms of increasing intra parcel correlation,
		// however it makes the most difference in the graph based analysis. 
		vnl_vector<float> g(x.size());
		for (int i = 0; i<x.size(); i++)
		{
		g[i] = dot_product(segmentSVDs[f-1],x[i]);
		}
		int indx = g.arg_max();
		segmentSVDs[f - 1] = x[indx];

		++f;
		itParcel.GoToBegin();
		itFMRI.GoToBegin();
	}

	endTime = std::clock();

	std::cout << "Finding Representative signals of parcels by SVD performed in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;

	std::cout << "Finding segment information done totally in " << double(endTime - startTimeTop) / CLOCKS_PER_SEC << " seconds" << std::endl;
}

Vav::Parcellation::Parcellation::ParcelImageType::Pointer Vav::Parcellation::Parcellation::GetNewParcelImage(){
	return finalParcelImage;
}

void Vav::Parcellation::Parcellation::findBoundaryVoxelsAndDisplacementField()
{
	/*
	This function finds 
	*/

	//    Performance calculation variables declared
	std::clock_t startTime, endTime;
	startTime = std::clock();
	std::cout << "Starting boundary voxel search and displacement field generation...\n";

	//    Temprorary variables
	DisplacementFieldType::PixelType dispVector;
	int oldLabelChecked = 0;
	int mainLabel = 0;
	float oldCorrelationValue = 0;
	float newCorrelationValue = 0;
	float cost = 0;
	int newLabel = 0;
	float deltaVolumeChange = 0;
    float volumeCost = 0;
    float correlationCost = 0;


	//    Declaration of displacement field region iterator
	typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementRegionIteratorType;
	DisplacementRegionIteratorType itDisplacement(displacementField, displacementField->GetLargestPossibleRegion());
	itDisplacement.GoToBegin();

	//    Displacement field reset
	dispVector.Fill(0);
	while (!itDisplacement.IsAtEnd()) {
		itDisplacement.Set(dispVector);
		++itDisplacement;
	}
	itDisplacement.GoToBegin();

	//    Declaration of neighborhood iterator
	typedef itk::ConstNeighborhoodIterator<ParcelImageType> ParcelImageConstNeighborhoodIteratorType;
	ParcelImageConstNeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	ParcelImageConstNeighborhoodIteratorType itNeighborParcel(radius, parcelImage, parcelImage->GetLargestPossibleRegion());
	itNeighborParcel.GoToBegin();

	//    Offsets
	std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType> offsets = { { { -1, 0, 0 } }, { { 1, 0, 0 } }, { { 0, -1, 0 } }, { { 0, 1, 0 } }, { { 0, 0, -1 } }, { { 0, 0, 1 } } };
	std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType>::iterator itOffsets;

	typedef itk::ImageRegionIterator<FMRIImageType> FMRIRegionIteratorType;
	FMRIRegionIteratorType itFMRI(fmriImage, fmriImage->GetLargestPossibleRegion());
	itFMRI.GoToBegin();

	typedef itk::ImageRegionIterator<DirectionOfEvolutionImageType> DirectionImageRegionIteratorType;
	DirectionImageRegionIteratorType itDirection(directionImage, directionImage->GetLargestPossibleRegion());
	itDirection.GoToBegin();

	while (!itNeighborParcel.IsAtEnd()) 
	{
		mainLabel = itNeighborParcel.GetCenterPixel();
		oldLabelChecked = 0;
		oldCorrelationValue = 0;
		newLabel = 0;
		dispVector.Fill(0);
		cost = 0;
		deltaVolumeChange = 0;
        volumeCost = 0;
        correlationCost = 0;

		if (mainLabel != 0) 
		{
			//oldCorrelationValue = dot_product(itFMRI.Get(), segmentSVDs[mainLabel - 1]);
			for (itOffsets = offsets.begin(); itOffsets != offsets.end(); itOffsets++) 
			{
				if (itNeighborParcel.GetPixel((*itOffsets)) != mainLabel && itNeighborParcel.GetPixel((*itOffsets)) != 0) 
				{
					newLabel = itNeighborParcel.GetPixel((*itOffsets));
					if (newLabel == oldLabelChecked) 
					{
						dispVector[0] = (*itOffsets)[0];
						dispVector[1] = (*itOffsets)[1];
						dispVector[2] = (*itOffsets)[2];
						itDisplacement.Value() += dispVector;
					}
					else 
					{
						newCorrelationValue = dot_product(itFMRI.Get(), segmentSVDs[newLabel - 1]);
						if (newCorrelationValue > oldCorrelationValue) {
							oldLabelChecked=newLabel;
							oldCorrelationValue = newCorrelationValue;
							dispVector[0] = (*itOffsets)[0];
							dispVector[1] = (*itOffsets)[1];
							dispVector[2] = (*itOffsets)[2];
							itDisplacement.Set(dispVector);
							correlationCost = dot_product(itFMRI.Get(), segmentSVDs[mainLabel - 1]) - newCorrelationValue;
						}
					}
				}
			}
			
			if (oldLabelChecked != 0) 
			{
                deltaVolumeChange = parcelVolumeChanges[mainLabel - 1] - parcelVolumeChanges[oldLabelChecked - 1];
                volumeCost = sigmoidFunction(deltaVolumeChange, volumeSigmoidConstant);
                volumeCost = -1 * (2 * volumeCost - 1);
                cost = (1 - lamda) * correlationCost + lamda * volumeCost;
            }

			itDisplacement.Value().Normalize();
			itDisplacement.Value().operator*=(cost * (-cost_coeff));
			itDirection.Set(cost);
		}
		else 
		{
			itDisplacement.Set(dispVector);
		}

		++itNeighborParcel;
		++itDisplacement;
		++itFMRI;
		++itDirection;
	}
	endTime = std::clock();
	std::cout << "Boundary voxel search and displacement field generation done in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
}



void Vav::Parcellation::Parcellation::smoothDisplacementField(){
	/*
	This function smooths the displacement field using a Gaussian kernel.
	*/
	std::cout << "Starting Gaussian smoothing of displacement field...\n";
	std::clock_t startTime, endTime;
	startTime = std::clock();
	typedef itk::DiscreteGaussianImageFilter<DisplacementFieldType, DisplacementFieldType> GaussianSmoothingFilterType;
	GaussianSmoothingFilterType::Pointer smoothingFilter = GaussianSmoothingFilterType::New();
	smoothingFilter->SetVariance(gaussianSmoothingFilterVariance);
	smoothingFilter->SetMaximumKernelWidth(gaussianSmoothingFilterMaxKernelSize);
	smoothingFilter->SetInput(displacementField);
	smoothingFilter->Update();
	displacementField = smoothingFilter->GetOutput();
	displacementField->Update();
	endTime = std::clock();
	std::cout << "Smoothing done in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
}


void Vav::Parcellation::Parcellation::repairDisplacementField(){
	/*
	Since we use the displacement field to warp the boundary image, we do not want jumps and foldings. 
	Therefore, at each iteration, boundary should not move more than one voxel so that it does not fold.
	In order to assure that, we set the limit for maximum norm of the vectors in the displacement field to one.
	*/
	std::cout << "Starting repairing displacement field...\n";
	std::clock_t startTime, endTime;
	startTime = std::clock();

	//    Declaration of displacement field region iterator
	typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementRegionIteratorType;
	DisplacementRegionIteratorType itDisplacement(displacementField, displacementField->GetLargestPossibleRegion());
	itDisplacement.GoToBegin();

	while (!itDisplacement.IsAtEnd()) {
		if (itDisplacement.Get().GetNorm() > 1) {
			itDisplacement.Value().Normalize();
		}
		++itDisplacement;
	}

	endTime = std::clock();
	std::cout << "Repairing done in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;

}


void Vav::Parcellation::Parcellation::warpParcelImage(){
	/*
	This function warps the boundary image (parcellation) using the displacement field generated by calculated cost,
	and nearest neighbor interpolation since we interpolate labels. During the warping, there might happen holes,
	so before proceeding, we check the interior part of the cortex if it has holes, if yes, we interpolate the holes. 
	We also check if the outside of the cortex is all zero. If not, we force it to be zero.
	*/

    std::cout << "Starting warping parcel image..\n";
    std::clock_t startTime, endTime;
    startTime = std::clock();

    typedef itk::WarpImageFilter<ParcelImageType, ParcelImageType, DisplacementFieldType> WarpingFilterType;
    WarpingFilterType::Pointer warpingFilter = WarpingFilterType::New();
    typedef itk::NearestNeighborInterpolateImageFunction<ParcelImageType,double> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
	warpingFilter->SetInput(parcelImage);
	warpingFilter->SetInterpolator(interpolator);
    warpingFilter->SetOutputSpacing(parcelImage->GetSpacing());
    warpingFilter->SetOutputOrigin(parcelImage->GetOrigin());
    warpingFilter->SetOutputDirection(parcelImage->GetDirection());
    warpingFilter->SetDisplacementField(displacementField);
    
    parcelImage = warpingFilter->GetOutput();
    parcelImage->Update();
    checkZeros();
//    Hole filling
	std::cout << numberOfHolesInTheCortex() << std::endl;

    while (numberOfHolesInTheCortex() > 0) 
    {
        growParcelRegionInTheCortex();
    }
    setOutsideOfCortexToZero();
    checkZeros();
    endTime = std::clock();

    std::cout << "Warping done in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
}


void Vav::Parcellation::Parcellation::checkZeros(){
	/*
	This function checks the interior and exterior part of the cortex to see 
	if cortex has holes and if non-cortex has non-zero labels.
	*/

    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();    

    ParcelImageRegionIteratorType itInitialParcel(initialParcelImage, initialParcelImage->GetLargestPossibleRegion());
    itInitialParcel.GoToBegin();
    
    int counter = 0;
    int toplamVoxel = 0;
    int disariTasan = 0;
    while (!itParcel.IsAtEnd()) 
    {
        if (itInitialParcel.Get() != 0) 
        {
            if (itParcel.Get() == 0) 
            {
                counter++;
            }
        }
        if (itParcel.Get() != 0) 
        {
            toplamVoxel++;
        }

        ++itParcel;
        ++itInitialParcel;
        if (itParcel.Get() != 0) 
        {
            if (itInitialParcel.Get() == 0)
             {
                disariTasan++;
            }
        }
    }

    for (int i = 0; i < 15; i++) 
    {
        std::cout << "=";
    }
    std::cout << std::endl;
    std::cout << "There are " << counter << " holes in cortex at this moment.\n";
    std::cout << "There are " << toplamVoxel << " voxels in the brain at this moment.\n";
    std::cout << "There are " << disariTasan << " voxels outside the cortex at this moment.\n";
    for (int i = 0; i < 15; i++) 
    {
        std::cout << "=";
    }
    std::cout << std::endl;
}


void Vav::Parcellation::Parcellation::setOutsideOfCortexToZero(){
	/*
	If the noncortex area has nonzero labels, this function makes them zero. These problems arise after warping the parcellation.
	*/

    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    ParcelImageRegionIteratorType itInitialParcel(initialParcelImage, initialParcelImage->GetLargestPossibleRegion());
    itInitialParcel.GoToBegin();

    while (!itParcel.IsAtEnd()) {

        if (itInitialParcel.Get() == 0) 
        {
            itParcel.Set(0);
        }
        ++itParcel;
        ++itInitialParcel;
    }
}

int Vav::Parcellation::Parcellation::numberOfHolesInTheCortex(){    

//    returns number of holes in the cortex using parcelImage and initialParcelImage

    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    ParcelImageRegionIteratorType itInitialParcel(initialParcelImage, initialParcelImage->GetLargestPossibleRegion());
    itInitialParcel.GoToBegin();
    int counter = 0;

    while (!itParcel.IsAtEnd()) 
    {
        if (itInitialParcel.Get() != 0) 
        {
            if (itParcel.Get() == 0) 
            {
                counter++;
            }
        }
        ++itParcel;
        ++itInitialParcel;
    }
    return counter;
}

void Vav::Parcellation::Parcellation::growParcelRegionInTheCortex(){

//    Temp parcel image declared

	ParcelImageType::Pointer tempParcel = ParcelImageType::New();

    ParcelImageType::RegionType::IndexType parcelStartIndex;
    parcelStartIndex[0] = 0;
    parcelStartIndex[1] = 0;
    parcelStartIndex[2] = 0;
    ParcelImageType::RegionType::SizeType parcelSize;
    parcelSize = parcelImage->GetLargestPossibleRegion().GetSize();
    ParcelImageType::RegionType parcelRegion;
    parcelRegion.SetIndex(parcelStartIndex);
    parcelRegion.SetSize(parcelSize);
    tempParcel->SetRegions(parcelRegion);
    tempParcel->Allocate();

//    Declaration of tempParcel region iterator

    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itTempParcel(tempParcel, tempParcel->GetLargestPossibleRegion());
    itTempParcel.GoToBegin();

    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();

    // parcelImage is copied
    while (!itParcel.IsAtEnd()){
        itTempParcel.Set(itParcel.Get());
        ++itTempParcel;
        ++itParcel;
    }
    itTempParcel.GoToBegin();
    itParcel.GoToBegin();

//    Declaration of initialParcelImage iterator
    ParcelImageRegionIteratorType itInitialParcel(initialParcelImage, initialParcelImage->GetLargestPossibleRegion());
    itInitialParcel.GoToBegin();

    //    Declaration of neighborhood iterator
    typedef itk::ConstNeighborhoodIterator<ParcelImageType> ParcelImageConstNeighborhoodIteratorType;
    ParcelImageConstNeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    ParcelImageConstNeighborhoodIteratorType itNeighborParcel(radius, parcelImage, parcelImage->GetLargestPossibleRegion());
    itNeighborParcel.GoToBegin();

    //    Offsets

    std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType> offsets = {{{-1,0,0}},{{1,0,0}},{{0,-1,0}},{{0,1,0}},{{0,0,-1}},{{0,0,1}}};
    std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType>::iterator itOffsets;

    unsigned int mainLabel = 0;

    while (!itInitialParcel.IsAtEnd()) 
    {
        if (itInitialParcel.Get() != 0) 
        {
//            inside the cortex now

            mainLabel = itNeighborParcel.GetCenterPixel();

            if (mainLabel == 0) 
            {
//                at a hole in the cortex
                for (itOffsets = offsets.begin(); itOffsets != offsets.end(); itOffsets++) 
                {
                    if (itNeighborParcel.GetPixel((*itOffsets)) != 0) 
                    {
                        itTempParcel.Set(itNeighborParcel.GetPixel((*itOffsets)));
                    }
                }
            }
        }
        ++itInitialParcel;
        ++itTempParcel;
        ++itNeighborParcel;
    }
    parcelImage = tempParcel;
    parcelImage->Update();
}

//========================================================================

//Volume things

void Vav::Parcellation::Parcellation::findInitialParcelVolumes(){
    initialNumberOfSegmentVoxels.set_size(numberOfSegments);
    initialNumberOfSegmentVoxels.fill(0);
    //    Culumating # of voxels and their fmri signals
    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    while (!itParcel.IsAtEnd()) {
        if (itParcel.Get() > 0) {
            initialNumberOfSegmentVoxels[itParcel.Get() - 1] += 1;
        }
        ++itParcel;
    }
    parcelVolumeChanges.clear();
    parcelVolumeChanges.set_size(numberOfSegments);
    parcelVolumeChanges.fill(0);
}

void Vav::Parcellation::Parcellation::findParcelVolumeChanges(){
    parcelVolumeChanges.clear();
    parcelVolumeChanges.set_size(numberOfSegments);
    parcelVolumeChanges.fill(0);
    for (unsigned int i = 0; i < numberOfSegments; i++) {
        parcelVolumeChanges[i] = (numberOfSegmentVoxels[i] - initialNumberOfSegmentVoxels[i]) / initialNumberOfSegmentVoxels[i];
    }
}

float Vav::Parcellation::Parcellation::sigmoidFunction(float x, float sigmoidConstant)
{
    return (1 / (1 + exp(-1 * sigmoidConstant * x)));
}

void Vav::Parcellation::Parcellation::SetResultsFilename(std::string a){
    resultFilename = a;
}
/*
void Vav::Parcellation::Parcellation::applyVolumeConstraintToDisplacementField(){
    
//    Declaration of displacement field region iterator
    typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementRegionIteratorType;
    DisplacementRegionIteratorType itDisplacement(displacementField, displacementField->GetLargestPossibleRegion());
    itDisplacement.GoToBegin();
    
    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    
    typedef itk::ImageRegionIterator<DirectionOfEvolutionImageType> DirectionImageRegionIteratorType;
    DirectionImageRegionIteratorType itDirection(directionImage, directionImage->GetLargestPossibleRegion());
    itDirection.GoToBegin();
    
    int segmentIndex = 0;
    float weight = 0;
    int counter = 0;
    
    while (!itParcel.IsAtEnd()) {
        
        if (itDirection.Get() != 0) {
            segmentIndex = itParcel.Get() - 1;
            if (segmentIndex < 0) {
                std::cout << "There is a problem with segment labels in the cortex!!!!\n";
                counter++;
            }
            weight = sigmoidFunction(parcelVolumeChanges[segmentIndex], volumeSigmoidConstant);
            if (weight > 0.5) {
                weight = 1 - weight;
            }
            if ((parcelVolumeChanges[segmentIndex] * itDirection.Get()) < 0) {
                weight = 1 - weight;
            }
            itDisplacement.Value() *= weight;
        }
        ++itParcel;
        ++itDirection;
        ++itDisplacement;
    }
//    std::cout << "Number of problems is " << counter << std::endl;
    
}
 */

//========================================================================

void Vav::Parcellation::Parcellation::calculate_and_switch()
{
    unsigned int oldLabelChecked = 0;
    unsigned int mainLabel = 0;
    float oldCorrelationValue = 0;
    float newCorrelationValue = 0;
    float cost = 0;
    unsigned int newLabel = 0;
    float benefit;
    
    typedef itk::ConstNeighborhoodIterator<ParcelImageType> ParcelImageConstNeighborhoodIteratorType;
    ParcelImageConstNeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    ParcelImageConstNeighborhoodIteratorType itNeighborParcel(radius, parcelImage, parcelImage->GetLargestPossibleRegion());
    itNeighborParcel.GoToBegin();
//    std::vector<float,>
    //    Offsets
    std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType> offsets = {{{-1,0,0}},{{1,0,0}},{{0,-1,0}},{{0,1,0}},{{0,0,-1}},{{0,0,1}}};
    std::vector<ParcelImageConstNeighborhoodIteratorType::OffsetType>::iterator itOffsets;
    
    typedef itk::ImageRegionIterator<FMRIImageType> FMRIRegionIteratorType;
    FMRIRegionIteratorType itFMRI(fmriImage, fmriImage->GetLargestPossibleRegion());
    itFMRI.GoToBegin();
    
    
    Vec2D cost_info;
    std::vector<float> row;
    cost_info.clear();

    while (!itNeighborParcel.IsAtEnd())
    {
        mainLabel = itNeighborParcel.GetCenterPixel();
        ParcelImageType::IndexType main_index = itNeighborParcel.GetIndex();
        
        newLabel = 0;
        benefit = 0;
        row.clear();
        
        if (mainLabel != 0)
        {
            oldCorrelationValue = dot_product(itFMRI.Get(), segmentSVDs[mainLabel - 1]);
            for (itOffsets = offsets.begin(); itOffsets != offsets.end(); itOffsets++)
            {
                if (itNeighborParcel.GetPixel((*itOffsets)) != mainLabel && itNeighborParcel.GetPixel((*itOffsets)) != 0)
                {
                    newLabel = itNeighborParcel.GetPixel((*itOffsets));
                    //new_index = itNeighborParcel.GetIndex(itOffsets);
                    newCorrelationValue = dot_product(itFMRI.Get(), segmentSVDs[newLabel - 1]);
                    
                    if (newCorrelationValue > oldCorrelationValue)
                    {   
                        row.clear();
                        benefit = newCorrelationValue - dot_product(itFMRI.Get(), segmentSVDs[mainLabel - 1]);
                        oldCorrelationValue = newCorrelationValue;
                        row.push_back(benefit);
                        row.push_back(newLabel);
                        row.push_back(main_index[0]);
                        row.push_back(main_index[1]);
                        row.push_back(main_index[2]);
                    }
                }
            }
            
            if(benefit>0)
            {
                cost_info.push_back(row);
            }
            
        }
        
        ++itNeighborParcel;
        ++itFMRI;
    }
    
    std::vector<int> indices(cost_info.size());
    for (unsigned int i = 0; i < cost_info.size(); i++){ indices[i] = i;}
    sort(indices.begin(), indices.end(), sort_indices(&cost_info));
	if (cost_info.size() != 0)
	{
		for (unsigned int i = 0; i < cost_info.size()/2; i++)
		{
			ParcelImageType::IndexType idx;
			idx[0] = cost_info[i][2];
			idx[1] = cost_info[i][3];
			idx[2] = cost_info[i][4];
			parcelImage->SetPixel(idx, cost_info[i][1]);

		}
	}
}


void Vav::Parcellation::Parcellation::fillEvaluationData(){

	std::cout << "Starting filling evaluation data...\n";
	std::clock_t startTime, endTime;
	startTime = std::clock();

	evaluationData.clear();

	for (unsigned int i = 0; i < numberOfSegments; i++) {
		unsigned int id = 0;
		itk::PointSet<vnl_vector<float>>::Pointer pointSet = itk::PointSet<vnl_vector<float>>::New();

		typedef itk::ImageRegionIterator<ParcelImageType> IteratorType;
		IteratorType it(parcelImage, parcelImage->GetLargestPossibleRegion());
		it.GoToBegin();

		typedef itk::ImageRegionIterator<FMRIImageType> FMRIIteratorType;
		FMRIIteratorType itFMRI(fmriImage, fmriImage->GetLargestPossibleRegion());
		itFMRI.GoToBegin();

		while (!it.IsAtEnd()) {
			if (it.Get() - 1 == i) {
				pointSet->SetPointData(id, itFMRI.Get());
				id++;
			}
			++it;
			++itFMRI;
		}

		evaluationData.push_back(pointSet);
	}

	endTime = std::clock();
	std::cout << "Evaluation data filled in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
}



void Vav::Parcellation::Parcellation::calculateMeanCorrelations(){

	std::cout << "Starting finding results...\n";
	std::clock_t startTime, endTime;
	startTime = std::clock();

	std::vector<float> temp;
	std::vector<int> temp_volume;

	for (unsigned int i = 0; i < numberOfSegments; i++) {
		//        std::cout << "Now, finding for segment " << i << std::endl;

		typedef itk::PointSet<vnl_vector<float>>::PointDataContainer PointDataContainer;
		PointDataContainer::Pointer pointData = PointDataContainer::New();
		pointData = evaluationData[i]->GetPointData();

		PointDataContainer::Iterator it1;
		PointDataContainer::Iterator it2;

		it1 = pointData->Begin();

		float sum = 0;
		double counter = 0;

		//        for (it1 = pointData->Begin(); it1 != pointData->End(); it1++) {
		//            for (it2 = it1; it2 != pointData->End(); it2++) {
		//                sum += double(dot_product(it1.Value(), it2.Value()));
		//                counter += 1;
		//                std::cout << counter << " ";
		//            }
		//        }

		for (it1 = pointData->Begin(); it1 != pointData->End(); it1++) {
			sum += float(dot_product(segmentMeans[i], it1.Value()));
			counter += 1;
		}

		sum /= counter;
		temp.push_back(sum);
		temp_volume.push_back(counter);

	}

	results.push_back(temp);
	results_volume.push_back(temp_volume);

	endTime = std::clock();
	std::cout << "Results found in " << double(endTime - startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;
}


//=================================================================
//Start of downsizing parcels

int downsizeSpecificParcel(itk::Image<unsigned int,3>::Pointer parcelImage, unsigned int classLabel, unsigned int newClassLabel){
    
    typedef itk::Image<unsigned int, 3> ParcelImageType;
    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    double counter = 0;
    
    typedef itk::Vector< double, 3 > MeasurementVectorType;
    typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
    SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize(3);
    MeasurementVectorType mv;
    vnl_vector<double> center(3);
    center.fill(0);
    while (!itParcel.IsAtEnd()) {
        if (itParcel.Get() == classLabel) {
            mv[0] = double(itParcel.GetIndex()[0]);
            mv[1] = double(itParcel.GetIndex()[1]);
            mv[2] = double(itParcel.GetIndex()[2]);
            sample->PushBack(mv);
            center.operator+=(mv.GetVnlVector());
            counter++;
        }
        ++itParcel;
    }
    
    center.operator*=(1/counter);
//    Printing something
//    std::cout << "There are " << counter << " voxels in parcel " << classLabel << std::endl;
//    std::cout << "Center is " << center << std::endl;
    
    typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample( sample );
    treeGenerator->SetBucketSize( 16 );
    treeGenerator->Update();
    
    typedef TreeGeneratorType::KdTreeType TreeType;
    typedef itk::Statistics::KdTreeBasedKmeansEstimator< TreeType > EstimatorType;
    EstimatorType::Pointer estimator = EstimatorType::New();
    EstimatorType::ParametersType initialMeans(6);
    initialMeans[0] = center[0];
    initialMeans[1] = center[1];
    initialMeans[2] = center[2];
    initialMeans[3] = center[0] + 0.1;
    initialMeans[4] = center[1] + 0.1;
    initialMeans[5] = center[2] + 0.1;
    estimator->SetParameters( initialMeans );
    estimator->SetKdTree( treeGenerator->GetOutput() );
    estimator->SetMaximumIteration( 200 );
    estimator->SetCentroidPositionChangesThreshold(0.0);
    estimator->StartOptimization();
    
    EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
    
    //    Printing estimated means
    /*
    for (unsigned int i = 0; i<2; i++) {
        std::cout << "Center " << i << " is ";
        for (int j = 0; j<3; j++) {
            std::cout << estimatedMeans[3*i + j] << "  ";
        }
        std::cout << std::endl;
    }
    */
    
    typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
    typedef itk::Statistics::MinimumDecisionRule DecisionRuleType; DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
    typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
    ClassifierType::Pointer classifier = ClassifierType::New();
    classifier->SetDecisionRule( decisionRule );
    classifier->SetInput( sample );
    classifier->SetNumberOfClasses( 2 );
    typedef ClassifierType::ClassLabelVectorObjectType ClassLabelVectorObjectType;
    typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
    typedef ClassifierType::ClassLabelType ClassLabelType;
    ClassLabelVectorObjectType::Pointer classLabelsObject = ClassLabelVectorObjectType::New();
    ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();
    ClassLabelType class1 = 1;
    classLabelsVector.push_back( class1 );
    ClassLabelType class2 = 2;
    classLabelsVector.push_back( class2 );
    classifier->SetClassLabels( classLabelsObject );
    
    typedef ClassifierType::MembershipFunctionVectorObjectType MembershipFunctionVectorObjectType;
    typedef ClassifierType::MembershipFunctionVectorType MembershipFunctionVectorType;
    MembershipFunctionVectorObjectType::Pointer membershipFunctionVectorObject = MembershipFunctionVectorObjectType::New();
    MembershipFunctionVectorType& membershipFunctionVector = membershipFunctionVectorObject->Get();
    int index = 0;
    for (unsigned int i = 0; i < 2; i++){
        MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
        MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
        for ( unsigned int j = 0; j < sample->GetMeasurementVectorSize(); j++ ){
            centroid[j] = estimatedMeans[index++];
        }
        membershipFunction->SetCentroid( centroid );
        membershipFunctionVector.push_back( membershipFunction.GetPointer() );
    }
    classifier->SetMembershipFunctions( membershipFunctionVectorObject );
    classifier->Update();
    
    
    ParcelImageType::IndexType indexOfParcelImage;
    int a = 0;
    int b = 0;
    const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
    ClassifierType::MembershipSampleType::ConstIterator iter = membershipSample->Begin();
    while ( iter != membershipSample->End() ) {
        if (iter.GetClassLabel() == 1) {
            a++;
        }
        if (iter.GetClassLabel() == 2) {
            b++;
            indexOfParcelImage[0] = (unsigned int) (iter.GetMeasurementVector()[0]);
            indexOfParcelImage[1] = (unsigned int) (iter.GetMeasurementVector()[1]);
            indexOfParcelImage[2] = (unsigned int) (iter.GetMeasurementVector()[2]);
            parcelImage->SetPixel(indexOfParcelImage, newClassLabel);
        }
        ++iter;
    }
//    std::cout << "Class a has " << a << " elements.\n";
//    std::cout << "Class b has " << b << " elements.\n";
    return a;
}

void Vav::Parcellation::Parcellation::SetMaxParcelPopulation(int maxNumber){
    maxNumberOfParcelPopulation = maxNumber;
}

void Vav::Parcellation::Parcellation::downsizeParcels(){
    
    std::clock_t startTime = std::clock();
    
    typedef itk::ImageRegionIterator<ParcelImageType> ParcelImageRegionIteratorType;
    ParcelImageRegionIteratorType itParcel(parcelImage, parcelImage->GetLargestPossibleRegion());
    itParcel.GoToBegin();
    
    labelsAfterDownsizing.clear();
    for (unsigned int i = 0; i < labels.size(); i++) {
        labelsAfterDownsizing.push_back(i+1);
    }
    unsigned int newClassLabel = labels.size() + 1;
    
    int index = 0;
    while (index < labelsAfterDownsizing.size()) {
        
//        Find number of voxels in parcel (index + 1)
        itParcel.GoToBegin();
        int numberOfVoxels = 0;
        while (!itParcel.IsAtEnd()) {
            if (itParcel.Get() == index + 1) {
                numberOfVoxels++;
            }
            ++itParcel;
        }
        
        while (numberOfVoxels > maxNumberOfParcelPopulation) {
            numberOfVoxels = downsizeSpecificParcel(parcelImage, index + 1, newClassLabel);
            newClassLabel++;
            labelsAfterDownsizing.push_back(index+1);
        }
        index++;
    }
    std::clock_t endTime = std::clock();
    numberOfSegments = labelsAfterDownsizing.size();
    std::cout << "Downsizing parcels with upper limit " << maxNumberOfParcelPopulation << " is done in "
    << double(endTime-startTime)/CLOCKS_PER_SEC << " seconds." << std::endl;
}


//End of downsizing parcels

//=================================================================


void Vav::Parcellation::Parcellation::doSomethingg()
{
	//downsizeParcels();
	findInitialParcelVolumes();
    findSegmentInformation();
    findBoundaryVoxelsAndDisplacementField();
    smoothDisplacementField();
	repairDisplacementField();
    fillEvaluationData();
    calculateMeanCorrelations();
    warpParcelImage();
    

	for (int i = 1; i < iteration_limit; i++)
	{
		
	findSegmentInformation();
    findParcelVolumeChanges();
    findBoundaryVoxelsAndDisplacementField();
    smoothDisplacementField();
	repairDisplacementField();
    fillEvaluationData();
    calculateMeanCorrelations();
    warpParcelImage();

		//calculate_and_switch();

		std::cout << "iteration : " << i << std::endl;
	}
	fillEvaluationData();
    calculateMeanCorrelations();
	createFinalParcelImage();
	std::cout << " Total Number of Parcels is " << numberOfSegments << std::endl;
	
	std::ofstream myFile;
    myFile.open(resultFilename);
    myFile << numberOfSegments << " " << iteration_limit << " ";
    for (int i = 0; i < results.size() ; i++) {
        for (int j = 0; j < results[i].size(); j++) {
            myFile << results[i][j] << " ";
        }
    }
    for (int i = 0; i < results_volume.size() ; i++) {
        for (int j = 0; j < results_volume[i].size(); j++) {
            myFile << results_volume[i][j] << " ";
        }
    }
    myFile.close();

}





