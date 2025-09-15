/******************************************************************************
 * @file VtkWriter.cpp
 * @brief Implementation of VTK output writer for CFD results
 *****************************************************************************/

#include "VtkWriter.hpp"

// VTK includes
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>

#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>

namespace VtkWriter
{

void writeVtkFile
(
    const std::string& filename,
    const std::vector<Vector>& allNodes,
    const std::vector<Face>& allFaces,
    const std::map<std::string,
    const ScalarField*>& scalarCellFields
) 
{
    // Create a vtkPolyData object
    vtkSmartPointer<vtkPolyData> polyData = 
        vtkSmartPointer<vtkPolyData>::New();
    
    // Create vtkPoints from the allNodes vector
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints(allNodes.size());
    
    // Add all nodes to the points
    for (size_t i = 0; i < allNodes.size(); ++i) 
    {
        points->SetPoint(i, allNodes[i].x(), allNodes[i].y(), allNodes[i].z());
    }
    
    // Set the points to the polyData
    polyData->SetPoints(points);
    
    // Create vtkCellArray to store the polygons (faces)
    vtkSmartPointer<vtkCellArray> polygons =
        vtkSmartPointer<vtkCellArray>::New();
    
    // Keep track of owner cell for each face for mapping cell data to face data
    std::vector<size_t> faceOwnerCells;
    faceOwnerCells.reserve(allFaces.size());
    
    // Add all faces as polygons
    for (const auto& face : allFaces) 
    {
        // Create a polygon for this face
        vtkSmartPointer<vtkPolygon> polygon = 
            vtkSmartPointer<vtkPolygon>::New();

        polygon->GetPointIds()->SetNumberOfIds(face.nodeIndices().size());
        
        // Set the node indices for this polygon
        for (size_t i = 0; i < face.nodeIndices().size(); ++i)
        {
            polygon->GetPointIds()->SetId(i, face.nodeIndices()[i]);
        }
        
        // Add the polygon to the cell array
        polygons->InsertNextCell(polygon);
        
        // Store the owner cell index for this face
        faceOwnerCells.push_back(face.ownerCell());
    }
    
    // Set the polygons to the polyData
    polyData->SetPolys(polygons);
    
    // Add scalar cell fields to the polyData
    // Since PolyData uses faces as cells, we need to map cell data to face data
    for (const auto& [fieldName, scalarField] : scalarCellFields)
    {
        if (!scalarField) continue;
        
        // Create a vtkDoubleArray for this field
        vtkSmartPointer<vtkDoubleArray> dataArray = 
            vtkSmartPointer<vtkDoubleArray>::New();

        dataArray->SetName(fieldName.c_str());
        dataArray->SetNumberOfComponents(1);
        dataArray->SetNumberOfTuples(allFaces.size());
        
        // Map cell data to face data using the owner cell index
        for (size_t faceIdx = 0; faceIdx < allFaces.size(); ++faceIdx)
        {
            size_t ownerCellIdx = faceOwnerCells[faceIdx];
            
            // Check if the owner cell index is valid
            if (ownerCellIdx < scalarField->size())
            {
                dataArray->SetValue(faceIdx, (*scalarField)[ownerCellIdx]);
            }
            else
            {
                // Set a default value if the owner cell index is out of bounds
                dataArray->SetValue(faceIdx, 0.0);
                std::cerr   << "Warning: Owner cell index " << ownerCellIdx
                            << " is out of bounds for field " << fieldName 
                            << std::endl;
            }
        }
        
        // Add the data array to the polyData's cell data
        polyData->GetCellData()->AddArray(dataArray);
    }
    
    // Write the polyData to file using vtkXMLPolyDataWriter
    // Note: The output file should use .vtp extension for VTK PolyData XML format
    // If .vtk extension is used, the file will still be written but may cause confusion
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    
    // Set the data mode to ASCII for better readability (optional)
    // writer->SetDataModeToAscii();
    
    // Set the data mode to Binary for better performance and smaller file size
    writer->SetDataModeToBinary();
    
    // Enable compression for smaller file size
    writer->SetCompressorTypeToZLib();
    
    // Write the file
    writer->Write();
    
    std::cout   << "VTK PolyData file written successfully: " 
                << filename << std::endl;

    std::cout   << "  - Number of points: " << allNodes.size() << std::endl;

    std::cout   << "  - Number of polygons: " << allFaces.size() << std::endl;

    std::cout   << "  - Number of scalar fields: " << scalarCellFields.size()
                << std::endl;
}

} // namespace VtkWriter
