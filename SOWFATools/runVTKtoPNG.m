pathVTK = 'W:\OpenFOAM\marcusbecker-2.4.0\simulationCases\2021_Paper_Marcus\3turb_yaw_lowTI\postProcessing\postProcessing\sliceDataInstantaneous';
pathOutput = 'C:\Users\marcusbecker\surfdrive\PhD_Surf\01_Research\01_FLORIDyn\02_Matlab\FLORIDynDDC\Simulations\2021_3T_Case_yaw_IandI\Data\HHFlowField';
nameVTK = 'U_slice_horizontal.vtk';
VTKtoPNG(pathVTK,nameVTK,pathOutput,[2,11])
pathOutput = 'C:\Users\marcusbecker\surfdrive\PhD_Surf\01_Research\01_FLORIDyn\02_Matlab\FLORIDynDDC\Simulations\2021_3T_Case_yaw_IandI\Data\HHFlowFieldAveraged';
VTKtoAvePNG(pathVTK,nameVTK,pathOutput,[3.5,9.5])