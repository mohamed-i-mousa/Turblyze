/******************************************************************************
 * @file BoundaryData.cpp
 * @brief Implementation of boundary condition data storage and retrieval
 *****************************************************************************/

#include "BoundaryData.h"
#include <stdexcept>

void BoundaryData::setFixedValue(Scalar s_val) 
{
    type_ = BCType::FIXED_VALUE;
    scalarValue_ = s_val; 
    valueType_ = BCValueType::SCALAR;

    vectorValue_ = Vector();
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedValue(const Vector& v_val) 
{
    type_ = BCType::FIXED_VALUE;
    vectorValue_ = v_val; 
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    gradientType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(Scalar s_grad) 
{
    type_ = BCType::FIXED_GRADIENT;
    scalarGradient_ = s_grad; 
    gradientType_ = BCValueType::SCALAR;

    vectorValue_ = Vector();
    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setFixedGradient(const Vector& v_grad) 
{
    type_ = BCType::FIXED_GRADIENT;
    vectorGradient_ = v_grad; 
    gradientType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    valueType_ = BCValueType::UNDEFINED;
}

void BoundaryData::setZeroGradient() 
{
    type_ = BCType::ZERO_GRADIENT;
    scalarGradient_ = S(0.0);
    vectorGradient_ = Vector(S(0.0), S(0.0), S(0.0));

    valueType_ = BCValueType::UNDEFINED;
    gradientType_ = BCValueType::UNDEFINED;
}


void BoundaryData::setNoSlip() 
{
    type_ = BCType::NO_SLIP;
    vectorValue_ = Vector(S(0.0), S(0.0), S(0.0));
    valueType_ = BCValueType::VECTOR;

    scalarValue_ = S(0.0);
    gradientType_ = BCValueType::UNDEFINED;
}

Scalar BoundaryData::fixedScalarValue() const 
{
    if (type_ == BCType::FIXED_VALUE && valueType_ == BCValueType::SCALAR)
    {
        return scalarValue_;
    }
    
    throw std::runtime_error
    (
        "Attempted to get fixed scalar value, but BC is not set to "
        "FIXED_VALUE with SCALAR type."
    );
}

const Vector& BoundaryData::fixedVectorValue() const 
{
    if (type_ == BCType::FIXED_VALUE && valueType_ == BCValueType::VECTOR)
    {
        return vectorValue_;
    }

    if (type_ == BCType::NO_SLIP && valueType_ == BCValueType::VECTOR) 
    {
        return vectorValue_;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed vector value, but BC is not set to "
        "FIXED_VALUE/NO_SLIP with VECTOR type."
    );
}

Scalar BoundaryData::fixedScalarGradient() const 
{
    if 
    (
        type_ == BCType::FIXED_GRADIENT 
     && gradientType_ == BCValueType::SCALAR
    ) 
    {
        return scalarGradient_;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed scalar gradient, but BC is not set to "
        "FIXED_GRADIENT with SCALAR type."
    );
}

const Vector& BoundaryData::fixedVectorGradient() const 
{
    if 
    (
        type_ == BCType::FIXED_GRADIENT 
     && gradientType_ == BCValueType::VECTOR
    ) 
    {
        return vectorGradient_;
    }

    throw std::runtime_error
    (
        "Attempted to get fixed vector gradient, but BC is not set to "
        "FIXED_GRADIENT with VECTOR type."
    );
}
