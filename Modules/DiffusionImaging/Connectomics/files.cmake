set(CPP_FILES

  # DataStructures Connectomics
  IODataStructures/mitkConnectomicsNetwork.cpp
  IODataStructures/mitkConnectomicsConstantsManager.cpp

  # Rendering Connectomics
  Rendering/mitkConnectomicsNetworkMapper3D.cpp
  Rendering/mitkConnectomicsRenderingSchemeProperty.cpp
  Rendering/mitkConnectomicsRenderingEdgeFilteringProperty.cpp
  Rendering/mitkConnectomicsRenderingNodeFilteringProperty.cpp
  Rendering/mitkConnectomicsRenderingEdgeColorParameterProperty.cpp
  Rendering/mitkConnectomicsRenderingEdgeRadiusParameterProperty.cpp
  Rendering/mitkConnectomicsRenderingNodeColorParameterProperty.cpp
  Rendering/mitkConnectomicsRenderingNodeRadiusParameterProperty.cpp
  Rendering/mitkConnectomicsRenderingNodeThresholdParameterProperty.cpp
  Rendering/mitkConnectomicsRenderingEdgeThresholdParameterProperty.cpp
  Rendering/mitkConnectomicsEnumerationSubclassesSerializer.cpp

  # Algorithms Connectomics
  Algorithms/mitkConnectomicsNetworkCreator.cpp
  Algorithms/mitkConnectomicsHistogramBase.cpp
  Algorithms/mitkConnectomicsDegreeHistogram.cpp
  Algorithms/mitkConnectomicsShortestPathHistogram.cpp
  Algorithms/mitkConnectomicsBetweennessHistogram.cpp
  Algorithms/mitkConnectomicsHistogramCache.cpp
  Algorithms/mitkConnectomicsSyntheticNetworkGenerator.cpp
  Algorithms/mitkConnectomicsSimulatedAnnealingPermutationBase.cpp
  Algorithms/mitkConnectomicsSimulatedAnnealingPermutationModularity.cpp
  Algorithms/mitkConnectomicsSimulatedAnnealingManager.cpp
  Algorithms/mitkConnectomicsSimulatedAnnealingCostFunctionBase.cpp
  Algorithms/mitkConnectomicsSimulatedAnnealingCostFunctionModularity.cpp
  Algorithms/mitkConnectomicsStatisticsCalculator.cpp
  Algorithms/mitkConnectomicsNetworkConverter.cpp
  Algorithms/mitkConnectomicsNetworkThresholder.cpp
  Algorithms/mitkFreeSurferParcellationTranslator.cpp
)

set(H_FILES
  # Rendering Connectomics
  Rendering/mitkConnectomicsNetworkMapper3D.h
  Rendering/mitkConnectomicsRenderingProperties.h
  Rendering/mitkConnectomicsRenderingSchemeProperty.h
  Rendering/mitkConnectomicsRenderingEdgeFilteringProperty.h
  Rendering/mitkConnectomicsRenderingNodeFilteringProperty.h
  Rendering/mitkConnectomicsRenderingEdgeColorParameterProperty.h
  Rendering/mitkConnectomicsRenderingEdgeRadiusParameterProperty.h
  Rendering/mitkConnectomicsRenderingNodeColorParameterProperty.h
  Rendering/mitkConnectomicsRenderingNodeRadiusParameterProperty.h
  Rendering/mitkConnectomicsRenderingNodeThresholdParameterProperty.h
  Rendering/mitkConnectomicsRenderingEdgeThresholdParameterProperty.h


  # Datastructures Connectomics
  IODataStructures/mitkConnectomicsConstantsManager.h
  IODataStructures/mitkConnectomicsNetwork.h
  IODataStructures/mitkConnectomicsNetworkProperties.h

  # Algorithms Connectomics
  Algorithms/mitkConnectomicsNetworkCreator.h
  Algorithms/mitkConnectomicsHistogramBase.h
  Algorithms/mitkConnectomicsDegreeHistogram.h
  Algorithms/mitkConnectomicsShortestPathHistogram.h
  Algorithms/mitkConnectomicsBetweennessHistogram.h
  Algorithms/mitkConnectomicsHistogramCache.h
  Algorithms/mitkConnectomicsSyntheticNetworkGenerator.h
  Algorithms/mitkConnectomicsSimulatedAnnealingPermutationBase.h
  Algorithms/mitkConnectomicsSimulatedAnnealingPermutationModularity.h
  Algorithms/mitkConnectomicsSimulatedAnnealingManager.h
  Algorithms/mitkConnectomicsSimulatedAnnealingCostFunctionBase.h
  Algorithms/mitkConnectomicsSimulatedAnnealingCostFunctionModularity.h
  Algorithms/itkConnectomicsNetworkToConnectivityMatrixImageFilter.h
  Algorithms/mitkConnectomicsStatisticsCalculator.h
  Algorithms/mitkConnectomicsNetworkConverter.h
  Algorithms/BrainParcellation/mitkCostFunctionBase.h
  Algorithms/BrainParcellation/mitkRandomParcellationGenerator.h
  Algorithms/BrainParcellation/mitkRegionVoxelCounter.h
  Algorithms/mitkFreeSurferParcellationTranslator.h
  Algorithms/mitkCorrelationCalculator.h
)


