// Umut Küçükaslan
// Mehmet Onay

#ifndef Parcellation_h
#define Parcellation_h
#define cost_coeff 10
#define iteration_limit 11
#define gaussianSmoothingFilterVariance 0.4
#define gaussianSmoothingFilterMaxKernelSize 3
#define maxParcelPopulation 800
#define lamda 0
#define volumeSigmoidConstant 10

//Standart Libraries
#include <stdio.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <Eigen/SVD>


//ITK Libraries
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVariableLengthVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkWarpImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkListSample.h"
#include "itkPointSet.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkSampleClassifierFilter.h"

//VTK Libraries

namespace Vav {
    namespace Parcellation {
        
        class sort_indices
        {
        // function to sort the items in a 2D vector (in the first algorithm voxel labels changed according to the following rule: 
        // boundary voxels found and we check how much they are correlated with the current and neighbor parcel. Relative cost function is calculated and sorted for all boundary voxels.
        // Voxels decreasing the relative cost function (i.e. increasing the correlation more for a possible label change) change their label   
        public:
            typedef std::vector<std::vector<float>> Vec2D;
            sort_indices(Vec2D* parr) : mparr(parr) {}
            bool operator()(int i, int j) const { return (*mparr)[i][0] > (*mparr)[j][0]; }
            
        private:
            Vec2D* mparr;
            
        };
        
        class Parcellation {
            
        public:
            typedef std::vector<std::vector<float>> Vec2D;
            typedef itk::Image<unsigned int, 3> ParcelImageType;
            typedef itk::Image<float, 4> RawFMRIImageType;
            typedef vnl_vector<float> TimeSignalType;
            typedef itk::Image<TimeSignalType, 3> FMRIImageType;

            Parcellation ();
            ~Parcellation();
            
            
            void SetParcelImage(ParcelImageType::Pointer p);
            void SetFMRIImage(RawFMRIImageType::Pointer p);
            void SetBoundOfCortexLabels(const unsigned int lower, const unsigned int upper);
            void Update();
            void SetMaxParcelPopulation(int maxNumber);
            void SetResultsFilename(std::string a);
			ParcelImageType::Pointer GetNewParcelImage();
            std::vector<ParcelImageType::PixelType> GetLabels();
            ParcelImageType::Pointer GetInitialParcelImage();
            FMRIImageType::Pointer GetVectorFMRIImage();
            void calculate_and_switch();
            void doSomethingg();
			float sigmoidFunction(float x, float sigmoidConstant);
            
            //        protected:
            
            
        private:
            //            Declaration of input parcel and fmri images
            ParcelImageType::Pointer rawParcelImage = ParcelImageType::New();
            RawFMRIImageType::Pointer rawFMRIImage = RawFMRIImageType::New();
            
            //            Declaration of input related specifications
            //            cortex voxel labels are between 11000 and 13000
            unsigned int lowerBoundOfCortexLabels;
            unsigned int upperBoundOfCortexLabels;
            
            //            Declaration of parcel and fmri images to be used
            ParcelImageType::Pointer parcelImage = ParcelImageType::New();
            ParcelImageType::Pointer initialParcelImage = ParcelImageType::New();
			ParcelImageType::Pointer finalParcelImage = ParcelImageType::New();
            FMRIImageType::Pointer fmriImage = FMRIImageType::New();
            
            //            Old labels
            std::vector<ParcelImageType::PixelType> labels;
            unsigned int numberOfSegments;
			void createFinalParcelImage();
            
            //            Segment information variables and functions
            vnl_vector<int>   numberOfSegmentVoxels;
            std::vector<FMRIImageType::PixelType> segmentMeans;
            std::vector<FMRIImageType::PixelType> segmentSVDs;
            void findSegmentInformation();
            vnl_vector<float>   initialNumberOfSegmentVoxels;
            void findInitialParcelVolumes();
            vnl_vector<float> parcelVolumeChanges;
            void findParcelVolumeChanges();
            void applyVolumeConstraintToDisplacementField();
            
            //            Boundary voxels variables and functions
            typedef itk::Vector<float,3> DisplacementVectorType;
            typedef itk::Image<DisplacementVectorType, 3> DisplacementFieldType;
            DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
            typedef itk::Image<float,3> DirectionOfEvolutionImageType;
            DirectionOfEvolutionImageType::Pointer directionImage = DirectionOfEvolutionImageType::New();
            void findBoundaryVoxelsAndDisplacementField();
            
            
            //            Smoothing displacement field
            void smoothDisplacementField();

			//			  Repair displacement field
			void repairDisplacementField();

            //            Warping parcel image
            void warpParcelImage();
            int numberOfHolesInTheCortex();
            void growParcelRegionInTheCortex();
            void setOutsideOfCortexToZero();
            void checkZeros();

            //            Calculate correlations
            std::vector<itk::PointSet<vnl_vector<float>>::Pointer> evaluationData;
            void fillEvaluationData();
            void calculateMeanCorrelations();
            std::vector<std::vector<float>> results;
            std::vector<std::vector<int>> results_volume;
            std::string resultFilename = "/Users/Mehmet/Desktop/VAVLab/network/AD/svd/unnamed.txt";

             // Parcel downsizing: in case we want to limit the max voxel population of a parcel. Parcel volume variability is usually high in the existing atlas.
             // Further examination of functional homogeneity is possible with this parameter. Although it increases the intra parcel correlation of the signals, it is a questionable procedure and not justified here.
            int maxNumberOfParcelPopulation;
            void downsizeParcels();
			std::vector<ParcelImageType::PixelType> labelsAfterDownsizing;
            
        };
        
    }
}

#endif /* Parcellation_hpp */
